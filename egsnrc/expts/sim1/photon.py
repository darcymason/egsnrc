import numpy

from egsnrc.vect import *  # XXX eventually do selective import
from egsnrc import random

def toy(particles, key):
    key, ran_floats = random.floats_0_1(key, len(particles.energy))
    energy = particles.energy - 5.0 * ran_floats
    status_mask = energy <= 0
    particles = particles.cut(status_mask)
    particles = particles.set(energy=energy[status_mask])

    return particles, key


def photon(particles, regions, media, scoring, howfar, ausgab):
    """Called from `shower` to track photons

    Parameters
    ----------
    particles: Particles
        Particles class with all particle state contained

    media: Media
        Media class with access to medium information

    howfar:  function
        Callback to `howfar` user-code method

    ausgab:  function
        Callback to `ausgab` user-code method for tracking events and energy

    Returns
    -------
    ircode: int
        1 for a normal return
    """

    PCUT_DISCARD = 1
    USER_PHOTON_DISCARD = 2
    PAIR_ELECTRONS_KILLED = 3

    while len(particles):  # one "step" for each photon, either interacting or reaching boundary
        # if eig <= pcut[irl_m1]: particle_outcome = PCUT_DISCARD
        pcut_mask = particles.energy < pcut[particles.region]
        ialive[STATUS, pcut_mask] = PCUT_DISCARD

        # if wt == 0.0:  particle_outcome = USER_PHOTON_DISCARD
        wt0_mask = falive[WT] == 0
        ialive[STATUS, wt0_mask] = USER_PHOTON_DISCARD

        # ASSUME (for now) that all particles at top of loop are alive
        #    later in loop should remove any that are not
        # alive_mask = iparticles[STATUS] == 0  # could use ialive[STATUS] if don't add particles
        # falive = fparticles[:, alive_mask]
        # ialive = iparticles[:, alive_mask]

        gle = numpy.log(falive[ENERGY])

        # Sample number of mfp to transport before interacting
        rnno35 = rng.random(falive.shape[1])
        rnno35[rnno35==0] = 1.0e-30  # could use `numpy.clip` but not quite same as old code
        dpmfp = -numpy.log(rnno35)

        # XXX do not special case vacuum (at least for now) - just use very low rho
        # if medium != 0:
        # non_vac = (ialive[MEDIUM] != 0)
        # non_vac_mediums = ialive[MEDIUM, non_vac]
        # # set interval gle, ge;
        # # lgle is integer index into upcoming arrays
        # lgle = (
        #     ge1[non_vac_mediums]*gle + ge0[non_vac_mediums]
        # ).astype(numpy.int32)

        media = ialive[MEDIUM]
        # Calc integer index into sigma arrays
        lgle = (ge1[media] * gle + ge0[media]).astype(numpy.int32)


        # evaluate gmfpr0 using gmfp(gle)
        # Appears to do (if not sin)
        # {P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM);

        # gmfpr0, 1 set in egs_scale_photon_xsection() in EGSnrc orig
        gmfpr0 = gmfp1[MEDIUM, lgle] * gle + gmfp0[MEDIUM, lgle]




def howfar(fparticles, iparticles):
    logger.debug("In howfar")

def ausgab(fparticles, iparticles):
    logger.debug("In ausgab")

if __name__ == "__main__":
    # Toy examples of arrays

    # Particle indexes
    PARTICLE_FLOAT_PARAMS = 11
    ENERGY, X, Y, Z, U, V, W, WT, DNEAR, TSTEP, USTEP = range(PARTICLE_FLOAT_PARAMS)

    PARTICLE_INT_PARAMS = 3
    STATUS, REGION, MEDIUM = range(PARTICLE_INT_PARAMS)

    NUM_PARTICLES = 5

    # float and int state of particles
    # Different energies, one with weight 0
    fphotons = numpy.array(
        [
            (20, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0),
            ( 2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0),
            (18, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),  # wt= 0
            (12, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),  # v=w = 1
        ],
        dtype=numpy.float64
    ).transpose()
    # photons in several regions, different mediums
    # could look up medium by region, but easier to just track it
    iphotons = numpy.array(
        [
            (0, 2, 1),
            (0, 1, 1),
            (0, 3, 2),
            (0, 4, 3),
        ],
        dtype=numpy.int32
    ).transpose()

    # Completely made-up data arrays for toy example
    pcut = numpy.array((99, 20, 15, 12, 2), dtype=numpy.float64) # by region
    ge1 = numpy.array((99, 1.1, 2.2, 3.3), dtype=numpy.float64) # by medium
    ge0 = numpy.array((99, 0.1, 0.2, 0.3), dtype=numpy.float64) # by medium

    print("\n\nStarting particles states")
    print(fphotons)
    print(iphotons)

    print("\nStarting 'physics data'")
    print("pcut")
    print(pcut)
    print("ge1, ge0")
    print(ge1)
    print(ge0)

    rng = numpy.random.default_rng(12345)

    photon(fphotons, iphotons, howfar, ausgab)