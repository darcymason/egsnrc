from egsnrc import egsfortran
from egsnrc.randoms import randomset
from egsnrc.params import *
from egsnrc.commons import *
from egsnrc.constants import *
from .photo import photo

from math import log

import logging
logger = logging.getLogger("egsnrc")


# EMPTY CALLBACKS ----
call_howfar_in_photon = None
electron_region_change = None
photon_region_change = None
photonuc_correction = None
photonuclear = None
rayleigh_correction = None
rayleigh_scattering = None
select_photon_mfp = None
start_new_particle = None

# The following callbacks allow the user to change the particle
# selection scheme (e.g., adding importance sampling (splitting,
# leading particle selection, etc.)).
# (default callback is `particle-selection-photon`
# which in turn has ``None`` replacement)
particle_selection_photon = None
particle_selection_compt = particle_selection_photon
particle_selection_pair = particle_selection_photon
particle_selection_photo = particle_selection_photon


def photon(howfar, ausgab):
    """Called from `shower` to track a photon

    Parameters
    ----------
    howfar:  function
        Callback to `howfar` user-code method

    ausgab:  function
        Callback to `ausgab` user-code method for tracking events and energy

    Returns
    -------
    ircode: int
        1 for a normal return
    """

    # $ comin_photon
    # COMIN/DEBUG,BOUNDS,MEDIA,MISC,EPCONT,PHOTIN,STACK,THRESH,
    #   UPHIOT,USEFUL,USER,RANDOM,EGS-VARIANCE-REDUCTION/

    # $ define_local_variables_photon
    # "Local PHOTON variables in order of their appearance"
    # $ENERGY PRECISION
    #     PEIG;   precise photon energy
    # ;
    # $REAL
    #     eig,    photon energy
    #     rnno35, random number for default MFP selection
    #     gmfpr0, photon MFP before density scaling and coherent correction
    #     gmfp,   photon MFP after density scaling
    #     cohfac, Rayleigh scattering correction
    #     rnno37, random number for Rayleigh scattering selection
    #     xxx,    random number for momentum transfer sampling in Rayleigh
    #     x2,     scaled momentum transfer in Rayleigh scattering event
    #     q2,     momentum transfer squared in Rayleigh scattering event
    #     csqthe, COSTHE**2
    #     rejf,   Rayleigh scattering rejection function
    #     rnnorj, random number for rejection in Rayleigh scattering
    #     rnno36, random number for interaction branching
    #     gbr1,   probability for pair production
    #     gbr2,   probability for pair + compton
    #     t,      used for particle exchange on the stack
    # Ali:photonuc, 2 lines
    #     photonucfac, photonuclear correction
    #     rnno39; random number for photonuclear selection (RNNO38 is taken)
    # ;
    # $INTEGER
    #     iarg,   parameter for AUSGAB
    #     idr,    parameter for AUSGAB
    #     irl,    region number
    #     lgle,   index for GMFP interpolation
    #     lxxx;   index for Rayleigh scattering cummulative distribution int.

    logger.debug("Entered PHOTON")
    ircode = 1 # set up normal return
    np_m1 = np - 1  # ** 0-based Python

    peig = e[np_m1]
    eig = peig  # energy of incident gamma

    irl = ir[np_m1]
    irl_m1 = irl - 1  # ** 0-based Python

    medium_m1 = medium - 1

    # Set up flags for discards (For Python, DLM 2021-02)
    particle_outcome = None
    PCUT_DISCARD = 1
    USER_PHOTON_DISCARD = 2
    PAIR_ELECTRONS_KILLED = 3


    # --- Inline replace: $ start_new_particle; -----
    if start_new_particle:
        start_new_particle()
    else:
        useful.medium = med[irl_m1]
    # End inline replace: $ start_new_particle; ----

    if eig <= pcut[irl_m1]:
        particle_outcome = PCUT_DISCARD

    # :PNEWENERGY:
    # Enter this loop for each photon with new energy
    while particle_outcome is None:
        logger.debug(":PNEWENERGY:")
        if wt[np_m1] == 0.0:  # added May 01
            particle_outcome = USER_PHOTON_DISCARD
            break # XXX

        gle = log(eig) # gle is gamma log energy

        #    here to sample no. mfp to transport before interacting

        # --- Inline replace: $ SELECT_PHOTON_MFP; -----
        if select_photon_mfp:
            select_photon_mfp()
        else:
            rnno35 = randomset()
            if rnno35 == 0.0:
                rnno35 = 1.e-30
            photin.dpmfp = -log(rnno35)
        # End inline replace: $ SELECT_PHOTON_MFP; ----
        # NOTE:  This template can also be over-ridden by other schemes,
        #        such as the 'exponential transform' technique.

        epcont.irold = ir[np_m1]  # Initialize previous region

        # :PNEWMEDIUM:
        while True:  # Here each time we change medium during photon transport
            logger.debug(":PNEWMEDIUM:")
            if medium != 0:
                # set interval gle, ge;
                lgle = int(ge1[medium_m1]*gle + ge0[medium_m1]) # Set pwlf interval
                lgle_m1 = lgle - 1  # ** 0-based
                # evaluate gmfpr0 using gmfp(gle)
                gmfpr0 = gmfp1[lgle_m1, medium_m1]*gle+ gmfp0[lgle_m1, medium_m1]

            # :PTRANS:
            INTERACTION_READY = False
            while True:  # photon transport loop
                if medium == 0:
                    epcont.tstep = vacdst
                else:
                    epcont.rhof = rhor[irl_m1] / rho[medium_m1]  # density ratio scaling template
                    gmfp = gmfpr0 / rhof
                    # --- Inline replace: $ RAYLEIGH_CORRECTION; -----
                    if rayleigh_correction:
                        rayleigh_correction()
                    else:
                        if iraylr[irl_m1] == 1:
                            # evaluate cohfac using cohe(gle)
                            cohfac = cohe1[lgle_m1, medium_m1]*gle+ cohe0[lgle_m1, medium_m1]
                            gmfp *= cohfac
                    # End inline replace: $ RAYLEIGH_CORRECTION; ----
                    # Ali:photonuc, 1 line
                    # --- Inline replace: $ PHOTONUC_CORRECTION; -----
                    if photonuc_correction:
                        photonuc_correction()
                    else:
                        if iphotonucr[irl_m1] == 1:
                            # evaluate photonucfac using photonuc(gle)
                            photonucfac == photonuc1[lgle_m1, medium_m1]*gle+ photonuc0[lgle_m1, medium_m1]
                            gmfp *= photonucfac
                    # End inline replace: $ PHOTONUC_CORRECTION; ----  # A PHOTONUCLEAR TEMPLATE
                    epcont.tstep = gmfp * dpmfp

                # Set default values for flags sent back from user
                epcont.irnew = ir[np_m1]  # set default new region number
                epcont.idisc = 0  # assume photon not discarded
                epcont.ustep = tstep  # transfer transport distance to user variable
                epcont.tustep = ustep

                # IF (USTEP > dnear[np_m1]) [;howfar();]
                # --- Inline replace: $ CALL_HOWFAR_IN_PHOTON; -----
                if call_howfar_in_photon:
                    call_howfar_in_photon()
                else:
                    if ustep > dnear[np_m1] or wt[np_m1] <= 0:
                        howfar()
                # End inline replace: $ CALL_HOWFAR_IN_PHOTON; ---- # The above is the default replacement

                # Now check for user discard request
                if idisc > 0:
                    # user requested immediate discard
                    particle_outcome = USER_PHOTON_DISCARD
                    break # XXX

                epcont.vstep = ustep # set variable for output code
                epcont.tvstep = vstep
                epcont.edep = pzero # no energy deposition on photon transport

                epcont.x_final = x[np_m1] + u[np_m1]*vstep
                epcont.y_final = y[np_m1] + v[np_m1]*vstep
                epcont.z_final = z[np_m1] + w[np_m1]*vstep

                if iausfl[TRANAUSB-1+1] != 0:
                    ausgab(TRANAUSB)

                # Transport the photon
                x[np_m1] = x_final
                y[np_m1] = y_final
                z[np_m1] = z_final
                # Deduct from distance to nearest boundary
                dnear[np_m1] -= ustep
                if medium != 0:
                    photin.dpmfp = max(0., dpmfp - ustep / gmfp) # deduct mfp's

                epcont.irold = ir[np_m1] # save previous region

                useful.medold = medium
                if irnew != irold:
                    # REGION CHANGE
                    # --- Inline replace: $ photon_region_change; -----
                    if photon_region_change:
                        photon_region_change()
                    else:
                        # --- Inline replace: $ electron_region_change; -----
                        if electron_region_change:
                            electron_region_change()
                        else:
                            ir[np_m1] = irnew
                            irl = irnew
                            irl_m1 = irl - 1
                            useful.medium = med[irl_m1]
                            medium_m1 - medium - 1
                        # End inline replace: $ electron_region_change; ----
                    # End inline replace: $ photon_region_change; ----

                # After transport call to user
                if iausfl[TRANAUSA-1+1] != 0:
                    ausgab(TRANAUSA)
                # oct 31 bug found by C Ma. PCUT discard now after AUSGAB call
                if eig <= pcut[irl_m1]:
                    particle_outcome = PCUT_DISCARD
                    break

                # Now check for deferred discard request. May have been set
                # by either howfar, or one of the transport ausgab calls
                if idisc < 0:
                    particle_outcome = USER_PHOTON_DISCARD
                    break

                if medium != medold:
                    break  # exit :PTRANS: loop

                if medium != 0 and dpmfp <= EPSGMFP:
                    # Time for an interaction
                    INTERACTION_READY = True
                    break  # EXIT :PNEWMEDIUM:

                # end :PTRANS: loop
            # --------------------------
            # in :PNEWMEDIUM: loop
            if INTERACTION_READY or particle_outcome is not None:
                break
            # end :PNEWMEDIUM: loop
        # ----------
        # in :PNEWENERGY: loop
        if particle_outcome is not None:
            break  # go to discard sections

        # continue PNEWENERGY loop with interaction ...

        # It is finally time to interact.
        # The following allows one to introduce rayleigh scattering
        # --- Inline replace: $ RAYLEIGH_SCATTERING; -----
        if rayleigh_scattering:
            rayleigh_scattering()
        else:
            if iraylr[irl_m1] == 1:
                rnno37 = randomset()
                if rnno37 <= (1.0 - cohfac):
                    if iausfl[RAYLAUSB-1+1] != 0:
                        ausgab(RAYLAUSB)
                    stack.npold = np
                    egsfortran.egs_rayleigh_sampling(
                        medium, e[np_m1], gle, lgle, costhe, sinthe
                    )
                    uphi(2,1)
                    if iausfl[RAYLAUSA-1+1] != 0:
                        ausgab(RAYLAUSA)
                    continue # goto :PNEWENERGY:
        # End inline replace: $ RAYLEIGH_SCATTERING; ----
        # Ali:photonuclear, 1 line
        # --- Inline replace: $ PHOTONUCLEAR; -----
        if photonuclear:
            photonuclear()
        else:
            if iphotonucr[irl_m1] == 1:
                rnno39 = randomset()
                if rnno39 <= (1.0 - PHOTONUCFAC):
                    if iausfl[PHOTONUCAUSB-1+1] != 0:
                        ausgab(PHOTONUCAUSB)
                    egsfortran.photonuc()
                    if iausfl[PHOTONUCAUSA-1+1] != 0:
                        ausgab(PHOTONUCAUSA)
                    continue  # :PNEWENERGY: loop
        # End inline replace: $ PHOTONUCLEAR; ----
        rnno36 = randomset() # This random number determines which interaction
        #    GBR1 = PAIR / (PAIR + COMPTON + PHOTO) = PAIR / GTOTAL

        # Evaluate gbr1 using gbr1(gle)
        gbr1 = gbr11[lgle_m1, medium_m1]*gle+ gbr10[lgle_m1, medium_m1]
        if rnno36 <= gbr1 and e[np_m1] > rmt2:
            # IT WAS A PAIR PRODUCTION
            if iausfl[PAIRAUSB-1+1] != 0:
                ausgab(PAIRAUSB)
            egsfortran.pair()
            if particle_selection_pair:
                particle_selection_pair()

            if iausfl[PAIRAUSA-1+1] != 0:
                ausgab(PAIRAUSA)

            if iq[np_m1] != 0:
                break  # EXIT :PNEWENERGY:
            else:  # this may happen if pair electrons killed via Russian Roul
                # :PAIR_ELECTRONS_KILLED:
                # If here, then gamma is lowest energy particle.
                peig = e[np_m1]
                eig = peig
                if eig < pcut[irl_m1]:
                    particle_outcome = PCUT_DISCARD
                    break
                continue  # repeat PNEWENERGY loop

        #     GBR2 = (PAIR + COMPTON) / GTOTAL
        # evaluate gbr2 using gbr2(gle)
        gbr2 = gbr21[lgle_m1, medium_m1]*gle+ gbr20[lgle_m1, medium_m1]
        if rnno36 < gbr2:
            # It was a compton
            if iausfl[COMPAUSB+1-1] != 0:
                ausgab(COMPAUSB)
            egsfortran.compt()
            if particle-selection-compt:
                particle-selection-compt()
            if iausfl[COMPAUSA+1-1] != 0:
                ausgab(COMPAUSA)
            if iq[np_m1] != 0:  # Not photon
                break  # XXX EXIT:PNEWENERGY:
        else:
            if iausfl[PHOTOAUSB+1-1] != 0:
                ausgab(PHOTOAUSB)
            photo()  # egsfortran.photo()
            if particle_selection_photo:
                particle_selection_photo()

            if np == 0 or np < npold:
                return ircode

            # The above may happen if Russian Roulette is on.
            # np < npold means that only electrons were created in the interaction
            # and that all of them were killed. Hence, the top particle on the
            # stack is from a previous interaction and may be in another region
            # To avoid problems with the :PNEWENERGY: loop logic, we simply force
            # a return to shower so that ELECTR or PHOTON are properly re-entered.
            # Changed by IK Dec. 21 2006 after D. Rogers and R. Taylor found a
            # wrong dose with brems splitting and Russian Roulette on in a low
            # energy calculation.

            if iausfl[PHOTOAUSA+1-1] != 0:
                ausgab(PHOTOAUSA)
            if iq[np_m1] != 0:
                break  # XXX EXIT :PNEWENERGY:
        # End of photo electric block

        # end :PNEWENERGY: LOOP ---------


    if particle_outcome is None:
        # If here, means electron to be transported next
        return ircode

    # ---------------------------------------------
    # Photon cutoff energy discard section
    # ---------------------------------------------
    if particle_outcome == PCUT_DISCARD:
        if medium > 0:
            if eig > ap[medium_m1]:
                idr = EGSCUTAUS
            else:
                idr = PEGSCUTAUS
        else:
            idr = EGSCUTAUS

        epcont.edep = peig  # get energy deposition for user
        # inline replace $ PHOTON-TRACK-END
        if iausfl[idr-1+1] != 0:
            ausgab(idr)
        # --- end inline replace
        ircode = 2
        stack.np -= 1
        return ircode

    # ---------------------------------------------
    # User requested photon discard section
    # ---------------------------------------------
    elif particle_outcome == USER_PHOTON_DISCARD:
        epcont.edep = peig
        if iausfl[USERDAUS-1+1] != 0:
            ausgab(USERDAUS)
        ircode = 2
        stack.np -= 1
        return ircode
    else:
        raise ValueError(f"Unhandled particle outcome ({particle_outcome})")
