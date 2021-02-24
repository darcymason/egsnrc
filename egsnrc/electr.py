from egsnrc.randoms import randomset
from egsnrc.params import *
from egsnrc.commons import *
from egsnrc.constants import *
from egsnrc.electr_steps import tstep_ustep
from .annih import annih
from .bhabha import bhabha
import logging
logger = logging.getLogger('egsnrc')

# EMPTY CALLBACKS ----
call_user_electron = None



electron_track_end = None

evaluate_bhabha_fraction = None
evaluate_ebrem_fraction = None
evaluate_pbrem_fraction = None

# The following allow the user to change the
# particle selection scheme (e.g., adding importance sampling
# such as splitting, leading particle selection, etc.).
# Default replacement is `particle_selection_electr` which in turn
# is None
particle_selection_electr = None
particle_selection_annih = particle_selection_electr
particle_selection_annihrest = particle_selection_electr
particle_selection_bhabha = particle_selection_electr
particle_selection_brems = particle_selection_electr
particle_selection_moller = particle_selection_electr

start_new_particle = None


ierust = 0 # To count negative ustep's

# Define ebrems as is used in several places
def ebrems(ausgab):
    if iausfl[BREMAUSB-1+1] != 0:  # ** 0-based
        ausgab(BREMAUSB)
    egsfortran.brems()

    if particle_selection_brems:
        particle_selection_brems()

    if iausfl[BREMAUSA-1+1] != 0:  # ** 0-based
        ausgab(BREMAUSA)


# ******************************************************************
#                                NATIONAL RESEARCH COUNCIL OF CANADA
def electr(hownear, howfar, ausgab) -> int:
    # ******************************************************************
    #    This subroutine has been almost completely recoded to include
    #    the EGSnrc enhancements.
    #
    #    Version 1.0   Iwan Kawrakow       Complete recoding
    #    Version 1.1   Iwan Kawrakow       Corrected implementation of
    #                                      fictitious method (important
    #                                      for low energy transport
    #    egsnrc v0.1   Darcy Mason  Transpiled to Python
    # ******************************************************************

    global ierust


    # --- Inline replace: $ CALL_USER_ELECTRON -----
    if call_user_electron:
        call_user_electron()
    # End inline replace: $ CALL_USER_ELECTRON ----

    np_m1 = np - 1  # ** 0-based

    ircode = 1
    """Set up normal return-which means there is a photon
    with less available energy than the lowest energy electron,
    so return to shower so it can call photon to follow it.
    (For efficiency's sake, we like to stay in this routine
        as long as there are electrons to process. That's why this
        apparently convoluted scheme of STACK contro is effected.)
    """

    epcont.irold = ir[np_m1]
    """Initialize previous region
    (ir[] is an integer that is attached to the particle's
    phase space. It contains the region
    number that the current particle is in.
    Np is the stack pointer, it points to where on the
    stack the current particle is.)
    """

    irl = irold  # region number in local variable
    irl_m1 = irl - 1

    # --- Inline replace: $ start_new_particle; -----
    if start_new_particle:
        start_new_particle()
    else:
        useful.medium = med[irl_m1]  # ** 0-based arrays
        medium_m1 = medium - 1
    # End inline replace: $ start_new_particle; ----


    # ************************************************************************
    while True:  # :NEWELECTRON: LOOP
        # Go once through this loop for each 'new' electron whose charge and
        # energy has not been checked

        # Save charge in local variable
        # (iq = -1 for electrons, 0 for photons and 1 for positrons)
        lelec = iq[np_m1]

        # XX below looks like float precision conversion, but latest egsnrc.macros
        # has both as real*8
        # "ENERGY PRECISION" is 'DOUBLE PRECISION' which appears to be real*8
        # eie is $REAL, also defined as real*8
        peie  = e[np_m1]  # precise energy of incident electron (double precision)
        eie   = peie  # energy incident electron (conversion to single)
        # eie   = numpy.array(peie, dtype=numpy.float32) # convert to single precision

        if eie <= ecut[irl_m1]:
            particle_outcome = ECUT_DISCARD
            break
            # (Ecut is the lower transport threshold.)

        # useful.medium = med[irl_m1] # (This renders the above assignment redundant!)
        # The above assignment is unnecessary, IK, June 2003

        if wt[np_m1] == 0.0:
            particle_outcome = USER_ELECTRON_DISCARD
            break

        # Follow the one particle  -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
        particle_outcome, lelke, eie, peie = tstep_ustep(
            lelec, medium, irl,
            eie, peie,
            hownear, howfar, ausgab)
        # -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
        if particle_outcome == LAMBDA_WARNING:
            ircode = lelke
            return ircode

        if particle_outcome in (USER_ELECTRON_DISCARD, ECUT_DISCARD):
            break  # from new-electron loop

        lelke_m1 = lelke - 1  # ** 0-based

        #  Now sample electron interaction
        if lelec < 0:
            # e-,check branching ratio
            # --- Inline replace: $ EVALUATE_EBREM_FRACTION; -----
            if evaluate_ebrem_fraction:
                ebr1 = evaluate_ebrem_fraction()
            else:
                # EVALUATE ebr1 USING ebr1(elke)
                ebr1 = ebr11[lelke_m1, medium_m1]*elke+ ebr10[lelke_m1, medium_m1]
            # End inline replace: $ EVALUATE_EBREM_FRACTION; ----

            rnno24 = randomset()
            # logger.info(f'random rnno24={rnno24}')
            if rnno24 <= ebr1:
                # It was bremsstrahlung
                ebrems(ausgab)
                np_m1 = np - 1 # if new particle
                if iq[np_m1] != 0:
                    continue  # new-electron loop
                return ircode  # Photon was selected, return to shower

            # It was Moller, but first check the kinematics.
            # However, if EII is on, we should still permit an interaction
            # even if E<moller threshold as EII interactions go down to
            # the ionization threshold which may be less than thmoll.
            if e[np_m1] <= thmoll[medium_m1] and eii_flag == 0:
                # (thmoll = lower Moller threshold)
                # Not enough energy for Moller, so
                # force it to be a bremsstrahlung---provided ok kinematically.
                if ebr1 <= 0: # Brems not allowed either.
                    continue  #  NEW-ELECTRON loop
                # It was bremsstrahlung
                ebrems(ausgab)
                np_m1 = np - 1 # if new particle
                if iq[np_m1] != 0:
                    continue  # new-electron loop
                return ircode  # Photon was selected, return to shower

            if iausfl[MOLLAUSB-1+1] != 0:  # ** 0-based
                ausgab(MOLLAUSB)
            egsfortran.moller()
            if particle_selection_moller:
                particle_selection_moller()

            if iausfl[MOLLAUSA-1+1] != 0:  # ** 0-based
                ausgab(MOLLAUSA)
            np_m1 = np - 1 # if new particle
            if iq[np_m1] == 0:
                return ircode

            continue  # NEW-ELECTRON loop Electron is lowest energy-follow it

        # e+ interaction. pbr1 = brems/(brems + bhabha + annih
        # --- Inline replace: $ EVALUATE_PBREM_FRACTION; -----
        if evaluate_pbrem_fraction:
            pbr1 = evaluate_pbrem_fraction()
        else:
            # EVALUATE pbr1 USING pbr1(elke)
            pbr1 = pbr11[lelke_m1, medium_m1]*elke+ pbr10[lelke_m1, medium_m1]
        # End inline replace: $ EVALUATE_PBREM_FRACTION; ----

        rnno25 = randomset()
        # logger.info(f'random rnno25={rnno25}')
        if rnno25 < pbr1:
            # It was bremsstrahlung
            ebrems(ausgab)
            np_m1 = np - 1 # if new particle
            if iq[np_m1] != 0:
                continue  # new-electron loop
            return ircode  # Photon was selected, return to shower

        # Decide between bhabha and annihilation
        # pbr2 is (brems + bhabha)/(brems + bhabha + annih)
        # --- Inline replace: $ EVALUATE_BHABHA_FRACTION; -----
        if evaluate_bhabha_fraction:
            prb2 = evaluate_bhabha_fraction()
        else:
            # EVALUATE pbr2 USING pbr2(elke)
            pbr2 = pbr21[lelke_m1, medium_m1]*elke+ pbr20[lelke_m1, medium_m1]
        # End inline replace: $ EVALUATE_BHABHA_FRACTION; ----

        if rnno25 < pbr2:
            # It is bhabha
            if iausfl[BHABAUSB-1+1] != 0:  # ** 0-based
                ausgab(BHABAUSB)
            bhabha()  # egsfortran.bhabha()
            if particle_selection_bhabha:
                particle_selection_bhabha()
            if iausfl[BHABAUSA-1+1] != 0:  # ** 0-based
                ausgab(BHABAUSA)
            np_m1 = np - 1 # if case of new particle
            if iq[np_m1] == 0:
                return ircode
        else:
            # It is in-flight annihilation
            if iausfl[ANNIHFAUSB-1+1] != 0:  # ** 0-based
                ausgab(ANNIHFAUSB)
            egsfortran.annih()
            if particle_selection_annih:
                particle_selection_annih()
            np_m1 = np - 1  # changing particles
            if iausfl[ANNIHFAUSA-1+1] != 0:  # ** 0-based
                ausgab(ANNIHFAUSA)

            return ircode  # EXIT NEW-ELECTRON loop, return to shower
            # After annihilation the gammas are bound to be the lowest energy
            # particles, so return and follow them.
        # end pbr2 else

    # loop on NEW-ELECTRON ***************************************************

    # Done NEW-ELECTRON loop, tally and return
    if particle_outcome not in (USER_ELECTRON_DISCARD, ECUT_DISCARD):
        return ircode
    # ---------------------------------------------
    # User requested electron discard section
    # ---------------------------------------------
    if particle_outcome == USER_ELECTRON_DISCARD:
        epcont.idisc = abs(idisc)

        if lelec < 0 or idisc == 99:
            epcont.edep =  e[np_m1] - prm
        else:
            epcont.edep =  e[np_m1] + prm

        if iausfl[USERDAUS-1+1] != 0:  # ** 0-based
            ausgab(USERDAUS)

        if idisc != 99:
            stack.np -= 1
            ircode = 2
            return ircode  # i.e., return to shower
        # else close up below
    # ---------------------------------------------
    # Electron cutoff energy discard section
    # ---------------------------------------------
    elif particle_outcome == ECUT_DISCARD:
        if medium > 0:
            if eie > ae[medium_m1]:
                idr = EGSCUTAUS
                if lelec < 0:
                    epcont.edep = e[np_m1] - prm
                else:
                    epcont.edep = peie - prm
            else:
                idr = PEGSCUTAUS
                epcont.edep =  e[np_m1] - prm
        else:
            idr = EGSCUTAUS
            epcont.edep =  e[np_m1] - prm

        # Define electron_track_end if you wish to modify the
        # treatment of track ends
        if electron_track_end:
            electron_track_end()
        else:
            iarg=idr
            if iausfl[iarg-1+1] != 0:  # ** 0-based
                ausgab(iarg)

    # in Mortran code, above two flags fall through into positron_annihilation
    #  and its following np and ircode setting

    if lelec > 0:  # positron_annihilation:  # NRCC extension 86/9/12
        # It's a positron. Produce annihilation gammas if edep < peie
        if edep < peie:
            if iausfl[ANNIHRAUSB-1+1] != 0:  # ** 0-based
                ausgab(ANNIHRAUSB)
            egsfortran.annih_at_rest()
            if particle_selection_annihrest:
                particle_selection_annihrest()

            if iausfl[ANNIHRAUSA-1+1] != 0:  # ** 0-based
                ausgab(ANNIHRAUSA)
            # Now discard the positron and take normal return to follow
            # the annihilation gammas.
            return ircode # i.e., return to shower

    # Final clean-up for any paths not exited some other way
    stack.np -= 1

    # tell shower an e- or un-annihilated e+ has been discarded
    ircode = 2

    return ircode  # i.e., return to shower
