from math import log, exp, sqrt, sin, cos
import sys
import numpy as np
from numba import cuda
import numba as nb

from numba import int32, float32
from numba.core.types import NamedTuple

import logging
from egsnrc.compton import compton

from egsnrc.config import KINT, device_jit, on_gpu
from egsnrc.constants import REST_MASS
from egsnrc.params import VACDST, EPSGMFP
from egsnrc.particles import replace_region_xyz
from egsnrc.particles import STEPPING, COMPTON, PHOTO
from egsnrc.particles import (
    INTERACTION_READY, GEOMETRY_DISCARD, PCUT_DISCARD, PHOTO_DISCARD, USER_PHOTON_DISCARD
)
from egsnrc.particles import set_particle
from egsnrc.media import GBR_PAIR, GBR_COMPTON
from egsnrc.photo import photo
from egsnrc import egsrandom

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)


@device_jit
def transport_photon(p, dpmfp, gle, regions, howfar):
    """Traverse regions as necessary according to the mfp

    Returns
    -------
    mod_p  Particle
        Copy of p but with new position and region, if applicable
    dpmfp   Float
        The remaining mfp
    status  int
        Flags for reason we have finished.
        E.g. INTERACTION_READY, or various ...DISCARD options
    """
    mod_p = p  # dummy in case of discards before mod_p created

    status = p.status
    # :PNEWMEDIUM:
    while status >= 0:  # Here each time we change medium during photon transport
        region= p.region
        medium = region.medium
        if medium.number != 0:
            # set interval gle, ge;
            lgle = KINT(medium.ge[1] * gle + medium.ge[0]) # Set pwlf interval
            # evaluate gmfpr0 using gmfp(gle)
            gmfpr0 = medium.gmfp01[1, lgle] * gle + medium.gmfp01[0, lgle]

        # :PTRANS:
        # Get region and medium index

        while True:  # photon transport loop
            if medium.number == 0:
                tstep = VACDST
            else:
                rhof = p.region.rho / medium.rho  # density ratio scaling template
                gmfp = gmfpr0 / rhof
                # XXX RAYLEIGH_CORRECTION; -----
                # XXX PHOTONUC_CORRECTION; -----
                tstep = gmfp * dpmfp

            # Set default values for flags sent back from user
            # tustep = ustep
            ustep, new_region, discard = howfar(p, regions, tstep)  # HOWFAR-----

            # Now check for user discard request
            if discard > 0:
                # user requested immediate discard
                status = USER_PHOTON_DISCARD
                break

            x = np.float32(p.x + p.u * ustep)
            y = np.float32(p.y + p.v * ustep)
            z = np.float32(p.z + p.w * ustep)

            # if iausfl[TRANAUSB-1+1] != 0:
            #     ausgab(TRANAUSB)

            # Transport the photon
            mod_p = replace_region_xyz(p, new_region, x, y, z)
            new_medium = new_region.medium
            # Deduct from distance to nearest boundary
            # dnear[np_m1] -= ustep
            if medium.number != 0:
                dpmfp = max(0., dpmfp - ustep / gmfp) # deduct mfp's

            if new_region.number != p.region.number:
                # REGION CHANGE
                medium = new_region.medium
                # XXX --- Inline replace: $ photon_region_change; -----

            # After transport call to user
            # if iausfl[TRANAUSA-1+1] != 0:
            #     ausgab(TRANAUSA)

            if mod_p.energy <= region.pcut:
                status = PCUT_DISCARD
                break

            # Now check for deferred discard request. May have been set
            # by either howfar, or one of the transport ausgab calls
            if discard < 0:
                status  = USER_PHOTON_DISCARD
                break

            if new_medium.number != p.region.medium.number:
                break  # exit :PTRANS: loop to medium loop

            if dpmfp <= EPSGMFP and new_medium.number != 0:
                # Time for an interaction
                status = INTERACTION_READY
                break  # EXIT :PTRANS: then :PNEWMEDIUM:

            p = mod_p
            # end :PTRANS: loop
            # --------------------------
        p = mod_p  # need to 'repeat' this in case break out of inner loop

    return mod_p, status, dpmfp, lgle

non_gpu_index = None

# Kernel
@cuda.jit
def photon_kernel(
    rng_states,
    num_particles, # XXX temp
    # iparticles, fparticles,
    regions, media, iscore, fscore
):
    """Main photon particle simulation kernel - implicitly called on each gpu thread

    Note
    ----
        The calling routine must set `egsnrc.photon.howfar` and `egsnrc.photon.ausgab`
        before calling this routine.
    """

    # Get unique grid index
    gid = cuda.grid(1) if on_gpu else non_gpu_index  # Global for testing only
    # if gid > len(fparticles):   # needed when have more GPU threads than particles
    if gid >= num_particles:
        return

    # CUDA Grid-Stride loop to reuse this thread
    threads_per_grid = cuda.blockDim.x * cuda.gridDim.x
    for i_arr in range(gid, num_particles, threads_per_grid):
        # Get the next particle from the source
        p = get_source_particle(rng_states, gid, regions)
        # p = set_particle(gid, regions, iparticles, fparticles)

        # Note Particle tuples are not mutable, so put `status` into mutable variable
        status = p.status

        if p.energy <= p.region.pcut:
            status = PCUT_DISCARD

        # :PNEWENERGY:
        # Follow the photon through all its interactions
        while status >= 0:  # Negative status for particle no longer tracked
            if p.w == 0.0:  # User photon discard
                status = USER_PHOTON_DISCARD
                break

            if p.energy <= p.region.pcut:
                status = PCUT_DISCARD
                break

            gle = log(p.energy)  # gamma log energy

            # Sample number of mfp to transport before interacting
            rnno35 = egsrandom.random_kfloat(rng_states, gid)
            if rnno35 == 0.0:
                rnno35 = 1.0e-30
            dpmfp = -log(rnno35)

            # print(f"{gle=} {rnno35=} {dpmfp=}")
            mod_p, status, dpmfp, lgle = transport_photon(
                p, dpmfp, gle, regions, howfar
            )

            ausgab(i_arr, status, p, mod_p, iscore, fscore) # if want to track change of position, region
            p = mod_p

            if status != INTERACTION_READY:
                break  # go to discard sections

            # It is finally time to interact.
            # The following allows one to introduce rayleigh scattering
            # XXX  --- Inline replace: $ RAYLEIGH_SCATTERING; -----
            # XXX --- Inline replace: $ PHOTONUCLEAR; -----

            # This random number determines which interaction
            rnno36 = egsrandom.random_kfloat(rng_states, gid)

            medium = mod_p.region.medium
            # Original Mortran
            # Evaluate gbr1 using gbr1(gle)
            #    GBR1 = PAIR / (PAIR + COMPTON + PHOTO) = PAIR / GTOTAL

            med_gbr = medium.gbr12[GBR_PAIR, :, lgle]
            gbr_pair = med_gbr[1] * gle + med_gbr[0]
            if rnno36 <= gbr_pair and p.energy > 2.0 * REST_MASS:
                # status = PAIR
                # mod_p = pair(rng_states, p)
                continue  # XXX for now skip pair
                # IT WAS A PAIR PRODUCTION
                # if iausfl[PAIRAUSB-1+1] != 0:
                #     ausgab(PAIRAUSB, ...)
                # p3, ...  = pair()
                # if particle_selection_pair:
                #     particle_selection_pair()

                # if iausfl[PAIRAUSA-1+1] != 0:
                #     ausgab(PAIRAUSA, ...)

                # if iq[np_m1] != 0:
                #     break  # EXIT :PNEWENERGY:
                # else:  # this may happen if pair electrons killed via Russian Roul
                #     # :PAIR_ELECTRONS_KILLED:
                #     # If here, then gamma is lowest energy particle.
                #     peig = e[np_m1]
                #     eig = peig
                #     if eig < pcut[irl_m1]:
                #         particle_outcome = PCUT_DISCARD
                #         break
                #     continue  # repeat PNEWENERGY loop

            #     GBR2 = (PAIR + COMPTON) / GTOTAL
            # evaluate gbr2 using gbr2(gle)
            med_gbr = medium.gbr12[GBR_COMPTON, :, lgle]
            gbr_compt = med_gbr[1] * gle + med_gbr[0]
            if rnno36 < gbr_compt:
                status = COMPTON
                mod_p = compton(rng_states, gid, p)
                # It was a compton
                # if iausfl[COMPAUSB+1-1] != 0:
                #     ausgab(COMPAUSB)
                # egsfortran.compt()
                # if particle-selection-compt:
                #     particle-selection-compt()
                # if iausfl[COMPAUSA+1-1] != 0:
                #     ausgab(COMPAUSA)
                # if iq[np_m1] != 0:  # Not photon
                #     break  # XXX EXIT:PNEWENERGY:
            else:
                # if iausfl[PHOTOAUSB+1-1] != 0: ...
                mod_p = photo(rng_states, gid, p)  # egsfortran.photo()
                status = PHOTO
                # if particle_selection_photo:...
                # if np == 0 or np < npold:
                #     return ircode
                # The above may happen if Russian Roulette is on....
                # if iausfl[PHOTOAUSA+1-1] != 0:
                #     ausgab(PHOTOAUSA)
            # End of photo electric block

            # Interaction done, looping back for next step
            ausgab(gid, status, p, mod_p, iscore, fscore)
            p = mod_p
            # end :PNEWENERGY: LOOP ---------

        # ---------------------------------------------
        # Photon cutoff energy discard section
        # ---------------------------------------------

        p = mod_p


        # if status == PCUT_DISCARD:
        #     if p.medium > 0:
        #         if eig > ap[medium_m1]:
        #             idr = EGSCUTAUS
        #         else:
        #             idr = PEGSCUTAUS
        #     else:
        #         idr = EGSCUTAUS

        #     epcont.edep = peig  # get energy deposition for user
        #     # inline replace $ PHOTON-TRACK-END
        #     if iausfl[idr-1+1] != 0:
        #         ausgab(idr)
        #     # --- end inline replace
        #     ircode = 2
        #     stack.np -= 1
        #     return ircode

        # ---------------------------------------------
        # User requested photon discard section
        # ---------------------------------------------
        # elif particle_outcome == USER_PHOTON_DISCARD:
        #     epcont.edep = peig
        #     if iausfl[USERDAUS-1+1] != 0:
        #         ausgab(USERDAUS)
        #     ircode = 2
        #     stack.np -= 1
        #     return ircode
        # else:
        #     raise ValueError(f"Unhandled particle outcome ({particle_outcome})")


def init(random_seed, num_particles):
    rng_states = create_xoroshiro128p_states(
        num_particles,
        seed=random_seed
    )
    return cuda.to_device(rng_states)


def run(particles, scoring_out, on_gpu=True):
    Py_major, Py_minor = sys.version_info.major, sys.version_info.minor
    print(f"Starting run with Numba {nb.__version__}, Python {Py_major}.{Py_minor}")
    print(f"Running {len(particles):,} particles")

    if cuda.is_available():
        print(cuda_details())
    else:
        thread_index = GridIterator()
        print("**** Running on CPU  ****")

    # To Debug on CPU:
    #  - set the flag True below
    # Add `gid` parameter to front of kernel call
    # Comment out the jit decorator for scoring function
    debug_on_cpu = False  # True
    from time import perf_counter


    particle_kernel.forall(len(fparticles))(
        dev_rng_states, dev_iparticles, dev_fparticles, dev_out
    )

    times = []
    for run in range(3):
        energies_np = np.random.random(num_photons) + 0.511 # catch both ko>2 and <2
        energies_np = energies_np.astype(np.float32)
        fparticles = np.zeros((num_photons, len(particle_fattrs)), dtype=np.float32)
        iparticles = np.zeros((num_photons, len(particle_iattrs)), dtype=np.int32)
        fparticles[:, ENERGY] = energies_np
        out = np.zeros(
            # 3 regions: 0, 1, and  2=escaped the geometry (just use eCompt)
            (num_photons, NUM_REGIONS + 1,  6),  # 6 for (nCompt, nPhoto, eCompt, ePhoto, nLost, eLost)
            dtype=np.float32
        )
        if not debug_on_cpu:
            dev_fparticles = cuda.to_device(fparticles)
            dev_iparticles = cuda.to_device(iparticles)
            dev_out = cuda.to_device(out)

        # Run
        if not debug_on_cpu:
            cuda.synchronize()
        start = perf_counter()
        if not debug_on_cpu:
            pass
            # Try to force typing
            # p = Particle(int32(0), int32(0), float32(1), float32(0))
            # mod_p = p
            # print(f"{type(iparticles[0, 0])=}")
            # print(f"{type(fparticles[0, 0])=}")

        else:
            egsrandom.random_kfloat = lambda r,i: np.random.random(1)
            for i in range(num_photons):
                particle_kernel.py_func(i, None, iparticles, fparticles, out)
        if not debug_on_cpu:
            cuda.synchronize()
        end = perf_counter()

        times.append(end - start)

    print("----------------------")
    print("Times:", ', '.join(f"{time_:>8.5} " for time_ in times), "seconds")
    if not debug_on_cpu:
        out = dev_out.copy_to_host()
    if len(out) < 50:
        print("Scoring - Particles/Region")
        print("nCompt     nPhoto     eCompt     ePhoto     nLost      eLost")
        print(out)
    for r in range(2):
        print("Energy Total by Region")
        print(f"Region {r}----")
        print("Compton: ", sum(out[:,r,SCORE_eCOMPTON]))
        print("Photo  : ", sum(out[:,r,SCORE_ePHOTO]))
    sum_compt = np.sum(out[:, :, SCORE_eCOMPTON])
    sum_photo = np.sum(out[:, :, SCORE_ePHOTO])
    sum_lost = np.sum(out[:, :, SCORE_eLost])
    print("Sum Compt, Photo, Lost", sum_compt, sum_photo, sum_lost)
    print("Energy in :", sum(fparticles[:, ENERGY]))
    print("Energy out:", sum((sum_compt, sum_photo, sum_lost)))

    # print("In particles: ---------------")
    # print(particles)
    # print("Out particles: ---------------")
    # print(out_particles)
    # sample = 100
    # print(f"First {sample} of various arrays:")
    # print(f"  Input energies: {energies_np[:sample]}")
    # print(f"  Input randoms : {rng_states_np[:sample]}")
    # print(f"  Energy out    : {host_energies[:sample]}")
    # print(f"  Costhe out    : {host_costhe[:sample]}")
