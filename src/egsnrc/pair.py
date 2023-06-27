from math import sqrt, log
from egsnrc import config
from egsnrc.constants import REST_MASS
import egsrandom

import logging

logger = logging.getLogger("egsnrc")


@config.device_jit
def set_e_angle(rng_states, gid, e_energy):
    pse = sqrt(max(0.0, (e_energy - REST_MASS) * (e_energy + REST_MASS)))
    costhe = egsrandom.random_kfloat(rng_states, gid)
    costhe = 1.0 - 2.0 * costhe
    sinthe = (
        REST_MASS
        * sqrt((1.0 - costhe) * (1.0 + costhe))
        / (pse * costhe + e_energy)
    )
    costhe = (e_energy * costhe + pse) / (pse * costhe + e_energy)

    return sinthe, costhe

@config.device_jit
def set_pair_angle(rng_states, gid, p):
    # if iprdst > 0:  -- here we always use 1
    for ichrg in (1, 2):
        if ichrg == 1:
            ese = energy_e1
        else:
            ese = ese2
        # if iprdst_use == 1: -- always 1
        if ichrg == 1:
            UPHI(2, 1)
        else:
            sinthe = -sinthe
            NP = NP + 1
            UPHI(3, 2)
    iq[np] = charge_e2
    iq[np - 1] = charge_e1
    return charge_e1, charge_e2, costhe, sinthe  # XXX ??, ......??


@config.device_jit
def pair(rng_states, gid, p):
    """Pair / Triplet production

    For a photon energy below 2.1 MeV, the energies of the pair
    particles are uniformly distributed in the allowed range.

    For a photon energy between 2.1 and 50 MeV the Bethe-Heitler
    cross section is employed, above 50 MeV the Coulomb-corrected
    Bethe-Heitler is used.

    Reference:  NRC report PIRS-701 (2021-04), sections 2.2.1.iv and 2.2.1.v

    Here we code only for the default IPRDST=1,
    leading order term for angular slection
    """

    medium = p.region.medium

    # if i_play_RR == 1:
    #     #  The user wants to play Russian Roulette. For pair
    #     #  it is much more efficient to do it BEFORE the
    #     #  actual sampling
    #     i_survived_RR = 0  # flag they all survive inititally
    #     if prob_RR <= 0:
    #         if n_RR_warning < MAX_RR_WARNING:
    #             n_RR_warning = n_RR_warning + 1
    #             logger.warning("Attempt to play Russian Roulette with prob_RR<0! ")
    #     else:
    #         rnno_RR = egsrandom.random_kfloat(rng_states, gid)
    #         if rnno_RR > prob_RR:  # The pair was killed
    #             i_survived_RR = 2  # flag both particles eliminated
    #             if np > 1:
    #                 np = np-1
    #             else:
    #                 #  get a proper exit from PHOTO, we have to leave at least
    #                 #  one particle on the stack
    #                 wt[np] = 0
    #                 e[np] = 0
    #             return
    #         else:
    #             wt[np] = wt[np] / prob_RR

    do_nrc_pair = False  # may be set True if itriplet > 0

    # if itriplet > 0 and p.energy > 4 * REST_MASS:
    #     itrip = dli_triplet * gle + bli_triplet
    #     ftrip = medium.a_triplet[itrip] * gle + medium.b_triplet[itrip]
    #     rnno34 = egsrandom.random_kfloat(rng_states, gid)
    #     if rnno34 < ftrip:
    #         #  Triplet production
    #         sample_triplet(rng_states, gid, p)
    #         return

    # if pair_nrc == 1:
    #     # Sample from the NRC pair cross section data base
    #     # (privided the energy is within the available range)
    #     k = p.energy / REST_MASS
    #     if k < nrcp_emax:
    #         do_nrc_pair = True
    #         if k <= nrcp_emin:
    #             ibin = 1;
    #         else:
    #             abin = 1 + log((k-2)/(nrcp_emin-2))*nrcp_dlei
    #             ibin = abin
    #             abin = abin - ibin
    #             rbin = egsrandom.random_kfloat(rng_states, gid)
    #             if rbin < abin:
    #                 ibin = ibin + 1

    #         xx = alias_sample1(NRC_PAIR_NX_1,nrcp_xdata,
    #                 nrcp_fdata(1,ibin,medium),nrcp_wdata(1,ibin,medium),
    #                 nrcp_idata(1,ibin,medium))
    #         #  The above returns the energy fraction of the positron
    #         if xx > 0.5:
    #             energy_e1 = prm*(1 + xx*(k-2))
    #             charge_e1 = 1
    #             energy_e2 = peig - energy_e1
    #             charge_e2 = -1
    #         else:
    #             energy_e2 = prm*(1 + xx*(k-2))
    #             charge_e2 = 1
    #             energy_e1 = peig - energy_e2
    #             charge_e1 = -1

    if not do_nrc_pair:
        # Calculate split of energy between the two particles
        if p.energy <= 2.1:
            # Below 2.1, use approximation
            # IK introduced this because uniform energy distribution
            # is probably a better approximation than a zero energy 'electron'
            # for low energy pair production

            rnno30 = egsrandom.random_kfloat(rng_states, gid)
            energy_e2 = REST_MASS + 0.5 * rnno30 * (p.energy - 2 * REST_MASS)
            energy_e1 = p.energy - energy_e2
        else:  # Above 2.1, must sample
            # Decide whether to use Bethe-Heitler or BH coulomb corrected
            if p.energy < 50.0:  # Use BH without Coulomb correction
                l = 5
                l1 = l + 1

                # Find the actual rejection maximum for this photon energy
                delta = 4 * medium.delcm / p.energy
                if delta < 1:
                    Amax = medium.dl1[l] + delta * (
                        medium.dl2[l] + delta * medium.dl3[l]
                    )
                    Bmax = medium.dl1[l1] + delta * (
                        medium.dl2[l1] + delta * medium.dl3[l1]
                    )
                else:
                    aux2 = log(delta + medium.dl6[l])
                    Amax = medium.dl4[l] + medium.dl5[l] * aux2
                    Bmax = medium.dl4[l1] + medium.dl5[l1] * aux2

                # and then calculate the probability for sampling from (br-1/2)**2
                aux1 = 1 - 2 * REST_MASS / p.energy
                aux1 = aux1 * aux1
                aux1 = aux1 * Amax / 3
                aux1 = aux1 / (Bmax + aux1)
            else:
                # Use BH Coulomb-corrected
                L = 7
                # The absolute maxima are close to the actual maxima at high energies
                # =>use the absolute maxima to save time
                # bpar: Prob. for the 12*(BR-1/2)**2 part in PAIR, eq. (2.7.105)
                Amax = medium.dl1[l]
                Bmax = medium.dl1[l + 1]
                aux1 = medium.bpar[2] * (1 - medium.bpar[1] * REST_MASS / p.energy)

            del0 = p.energy * medium.delcm
            Eavail = p.energy - 2 * REST_MASS

            while True:
                rnno30 = egsrandom.random_kfloat(rng_states, gid)
                rnno31 = egsrandom.random_kfloat(rng_states, gid)
                rnno34 = egsrandom.random_kfloat(rng_states, gid)
                if rnno30 > aux1:  # use the uniform part
                    br = 0.5 * rnno31
                    rejmax = Bmax
                    l1 = l + 1
                else:  # use the (br-1/2)**2 part of the distribution
                    rnno32 = egsrandom.random_kfloat(rng_states, gid)
                    rnno33 = egsrandom.random_kfloat(rng_states, gid)
                    br = 0.5 * (1 - max([rnno31, rnno32, rnno33]))
                    rejmax = Amax
                    l1 = l

                Eminus = br * Eavail + REST_MASS
                Eplus = p.energy - Eminus
                delta = del0 / (Eminus * Eplus)
                if delta < 1:
                    rejf = medium.dl1[l1] + delta * (
                        medium.dl2[l1] + delta * medium.dl3[l1]
                    )
                else:
                    rejf = medium.dl4[l1] + medium.dl5[l1] * log(delta + medium.dl6[l1])

                if rnno34 * rejmax <= rejf:
                    break

            energy_e2 = Eminus
            energy_e1 = Eplus
        # Have energies, now randomly determine charge
        rnno34 = egsrandom.random_kfloat(rng_states, gid)
        charge_e1 = -1 if rnno34 < 0.5 else 1
        charge_e2 = -charge_e1

    # Energy done, now calculated angles
    ESE2 = energy_e2
    e[np] = energy_e1
    E[NP + 1] = energy_e2
    #    This average angle of emission for both pair production and
    #    bremsstrahlung is much smaller than the average angle of
    #    multiple scattering for delta T transport=0.01 R.L.
    #    the initial and final momenta are coplanar
    #    set up a new 'electron'
    # --- Inline replace: $ SET_PAIR_ANGLE; -----
    charge_e1, charge_e2, costhe, sinthe  # XXX??, ......?? = set_pair_angle(eig, ...?)

    # End inline replace: $ SET_PAIR_ANGLE; ----

    #  DEFAULT FOR $ SET-PAIR-ANGLE; is to select the angle from the leading term
    #  of the angular distribution
    # UPHI(1,1)
    #    SET UP A NEW 'ELECTRON'
    # NP=NP+1
    # sinthe=-sinthe
    # UPHI(3,2)
    # iq[np]=charge_e2
    # IQ[NP-1]=charge_e1
