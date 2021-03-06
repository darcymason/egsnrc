from math import sqrt, exp
import sys

from .commons import *
from .params import *
from .angles import uphi
from .randoms import randomset
from .russian import russian_roulette

import logging
logger = logging.getLogger("egsnrc")

# EMPTY CALLBACKS ----
select_photoelectron_direction = None


def photo():
    """Photoelectic interaction

    Programmers:  I. Kawrakow, complete recoding,
                            Fluorescent X-rays, Auger,
                            Coster-Kronig treated in RELAX
                A.F. Bielajew (NRC) photoelectric angular distn

    2021-02 Darcy Mason - transpile to Python
    """

    logger.debug("Enter PHOTO")
    # comin_photo
    # COMIN/BOUNDS,DEBUG,EDGE,EGS-VARIANCE-REDUCTION,EPCONT,
    # MEDIA,PHOTIN,RANDOM,STACK,UPHIOT,USEFUL/

    # $ define_variables_for_select_photoelectron_direction
    # $REAL eelec,  total energy of photo-electron
    #       beta,   velocity of electron in units of c
    #       gamma,  total energy of photo-electron in units of rm
    #       alpha,  kinematic factor
    #       ratio,  =beta/alpha
    #       rnpht,  random number
    #       fkappa, aux. variable for costhe calculation
    #       xi,     used in rejection function calculation
    #       sinth2, sinthe**2

    # $ define_local_variables_photo
    # Local PHOTON variables in order of their appearance"
    # $ENERGY PRECISION
    #     peig;   "precise photon energy"
    # ;
    # $real
    #     eig,     photon energy
    #     rnno35,  random number for default mfp selection
    #     gmfpr0,  photon mfp before density scaling and coherent correction
    #     gmfp,    photon mfp after density scaling
    #     cohfac,  rayleigh scattering correction
    #     rnno37,  random number for rayleigh scattering selection
    #     xxx,     random number for momentum transfer sampling in rayleigh
    #     x2,      scaled momentum transfer in rayleigh scattering event
    #     q2,      momentum transfer squared in rayleigh scattering event
    #     csqthe,  costhe**2
    #     rejf,    rayleigh scattering rejection function
    #     rnnorj,  random number for rejection in rayleigh scattering
    #     rnno36,  random number for interaction branching
    #     gbr1,    probability for pair production
    #     gbr2,    probability for pair + compton
    #     t,       used for particle exchange on the stack
    #  Ali:photonuc, 2 lines
    #     photonucfac,  photonuclear correction
    #     rnno39;  random number for photonuclear selection (rnno38 is taken)
    # ;
    # $INTEGER
    #     iarg,    parameter for ausgab
    #     idr,     parameter for ausgab
    #     irl,     region number
    #     lgle,    index for gmfp interpolation
    #     lxxx;    index for Rayleigh scattering cummulative distribution int.


    n_warning = 0  # from Fortran Data statement

    if mcdf_pe_xsections:
        egsfortran.egs_shellwise_photo()
        return

    stack.npold = np # Set the old stack counter
    np_m1 = np - 1  # ** 0-based
    medium_m1 = medium - 1

    peig = e[np_m1]
    irl = ir[np_m1]
    irl_m1 = irl - 1  # ** 0-based Python

    if peig < edge_energies[2-1, 1-1]:  # ** 0-based for Python
        if n_warning < 100:
            n_warning += 1
            logger.info(
                f' Subroutine PHOTO called with E = {peig}'
                ' which is below the current min. energy of 1 keV! '
            )
            logger.info(' Converting now this photon to an electron, ')
            logger.info(' but you should check your code! ')

        iq[np_m1] = -1
        e[np_m1] = peig + prm
        return

    iZ = iedgfl[irl_m1]
    do_relax = False
    edep = pzero
    if iedgfl[irl_m1] != 0:
        #  User requested atomic relaxations
        #  first sample the element
        if nne[medium_m1] == 1:
            iZ = int( zelem[medium_m1, 1-1] + 0.5 )  # -1 for ** 0-based Python
            iZ_m1 = iZ - 1  # ** 0-based
            for j_m1 in range(edge_number[iZ_m1]):
                if peig >= edge_energies[j_m1, iZ_m1]:
                    break
            j = j_m1 + 1
        else:
            aux = peig*peig
            aux1 = aux*peig
            aux = aux*sqrt(peig)
            sigtot = 0
            for k_m1 in range(nne[medium_m1]):
                iZ = int( zelem[medium_m1, k_m1] + 0.5 )
                iZ_m1 = iZ - 1
                if iZ < 1 or iZ > MXELEMENT:
                    logger.info(' Error in PHOTO: ')
                    logger.critical(
                        '   Atomic number of element {k_m1 + 1}'
                    ' in medium {medium} is not between 1 and {MXELEMENT}'
                    )
                    sys.exit(-1)


                if peig > edge_energies(1-1, iZ_m1):

                    j = 1
                    sigma = (edge_a[1-1,iZ_m1] + edge_b[1-1,iZ_m1]/peig +
                        edge_c[1-1,iZ_m1]/aux + edge_d[1-1,iZ_m1]/aux1)/peig
                else:
                    for j_m1 in range(1, edge_number[iZ_m1]):
                        if peig >= edge_energies[j_m1, iZ_m1]:
                            j = j_m1 - 1  # do we need to update j?  XXX
                            break
                    j = j_m1 + 1

                    sigma = edge_a[j_m1, iZ_m1] + gle*(edge_b[j_m1, iZ_m1] + gle*(edge_c[j_m1, iZ_m1] +
                            gle*edge_d[j_m1, iZ_m1] ))
                    sigma = exp(sigma)

                sigma *= pz[medium_m1, k_m1]
                sigtot = sigtot + sigma
                probs[k_m1] = sigma
                ints[k_m1] = j

            br = randomset() * sigtot
            for k_m1 in range(nne[medium_m1]):
                br -= probs[k_m1]
                if br <= 0:
                    break

            iZ = int( zelem[medium_m1, k_m1] + 0.5 )
            iZ_m1 = iZ - 1
            j  = ints[k_m1]
            j_m1 = j - 1

        #  Now we know the atomic number (iZ) and the energy interval the
        #  photon energy is in (j). It is time to sample the shell the photon
        #  is interacting with.
        #  left for now as before, to be changed!!!
        if peig <= binding_energies[MXSHELL-1, iZ_m1]:  # ** -1 for 0-based
            # Outer shells, no atomic relaxation
            # EADL relax: Below  M2-shell -> just emit e-
            iq[np_m1] = -1
            e[np_m1] = peig + prm
        else:
            br = randomset()  # ftot = 1
            for k_m1 in range(MXINTER):
                if peig > binding_energies[k_m1, iZ_m1]:
                    if br < interaction_prob[k_m1, iZ_m1]:
                        break
                    br = (br - interaction_prob[k_m1, iZ_m1])/(1 - interaction_prob[k_m1, iZ_1])
            k = k_m1 + 1
            # Interaction possible with any shell from k=1 to MXSHELL
            # Defaults to MXSHELL interaction if for loop completes
            # ****************
            # EADL APPROACH 1: Do not allow interaction below L3. Deviates
            # **************** from previous EGSnrc approach as it doesn't
            #                  generate e- nor x-rays from <M> and <N> shells.
            if eadl_relax and k > 4:
                # No initial vacancy below L3 for now, just emit e-
                iq[np_m1] = -1
                e[np_m1] = peig + prm
            else:
                # EADL:    Interacts with K,L1..L3 shells
                # default: Interacts with K,L1..L3,<M>, and <N> shells
                e_vac = binding_energies[k_m1, iZ_m1]
                e[np_m1] = peig - e_vac + prm
                do_relax = True
                iq[np_m1] = -1
    else:
        e[np_m1] = peig + prm
        iq[np_m1] = -1

    if iq[np_m1] == -1:
        # --- Inline replace: $ SELECT_PHOTOELECTRON_DIRECTION; -----
        if select_photoelectron_direction:
            select_photoelectron_direction()
        else:
            #  ================================
            if iphter[ir[np_m1]] == 1:
                eelec = e[np_m1]
                if eelec > ecut[ir[np_m1]]:
                    beta = sqrt((eelec - rm) * (eelec + rm)) / eelec
                    gamma = eelec / rm
                    alpha = 0.5*gamma - 0.5 + 1./gamma
                    ratio = beta / alpha
                    while True:
                        rnpht = 2.*randomset() - 1.
                        if ratio <= 0.2:
                            fkappa = rnpht + 0.5*ratio*(1. - rnpht)*(1. + rnpht)
                            if gamma < 100:
                                uphiot.costhe = (beta + fkappa)/(1. + beta*fkappa)
                            else:
                                if fkappa > 0:
                                    uphiot.costhe = (
                                        1 - (1 - fkappa)*(gamma - 3)
                                        / (2*(1 + fkappa)*(gamma - 1)**3)
                                    )
                                else:
                                    uphiot.costhe = (
                                        (beta + fkappa)
                                        / (1. + beta*fkappa)
                                    )

                            # xi = 1./(1. - beta*costhe); <-- this numerically problematic
                            #                             at high energies, ik
                            xi = (1 + beta*fkappa)*gamma*gamma
                        else:
                            xi = gamma*gamma*(
                                1. + alpha*(
                                    sqrt(1. + ratio*(2.*rnpht + ratio)) - 1.
                                )
                            )
                            uphiot.costhe = (1. - 1./xi)/beta

                        sinth2 = max(0.,(1. - costhe)*(1. + costhe))

                        if randomset() <= 0.5*(1. + gamma)*sinth2*xi/gamma:
                            break  # exit loop

                    uphiot.sinthe = sqrt(sinth2)
                    uphi(2, 1)

        # End inline replace: $ SELECT_PHOTOELECTRON_DIRECTION; ---- # Samples photo-electron direction

    # ****************
    # EADL APPROACH 2: PE interactions with K, L1...L3,<M> and <N> shells,
    # **************** but vacancies below L3 deposit energy locally. It wont
    #                  produce x-rays from <M> and <N> shells.
    # IF ($ EADL_RELAX and k > 4)[
    #    edep = e_vac; do_relax = False
    # ]
    if do_relax:
        egsfortran.relax(e_vac, k, iZ)

    if edep > 0:
        if (iausfl[PHOTXAUS-1+1] != 0):
            ausgab(PHOTXAUS)  # generates IARG = 4 call


    russian_roulette(npold-1)  # **0-based  -1

