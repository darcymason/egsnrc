
from egsnrc.config import device_jit

import logging
logger = logging.getLogger("egsnrc")

# EMPTY CALLBACKS ----
select_low_energy_pair_prodiction = None
set_pair_angle = None
set_pair_rejection_function = None


@device_jit
def select_low_energy_pair_production(peig):
#   $RANDOMSET RNNO30; $RANDOMSET rnno34;
#   PESE2 = PRM + 0.5*RNNO30*(PEIG-2*PRM); PESE1 = PEIG - PESE2;
#   IF( rnno34 < 0.5 ) [ iq1 = -1; iq2 = 1; ] ELSE [ iq1 = 1; iq2 = -1; ]
# }
# " IK introduced this macro because uniform energy distribution"
# " is probably a better approximation than a zero energy 'electron'"
# " for low energy pair production"

    rnno30 = randomset()
    rnno34 = randomset()
    pese2 = REST_MASS + 0.5*rnno30*(peig-2*REST_MASS)
    pese1 = peig - pese2
    if rnno34 < 0.5:
        iq1 = -1
        iq2 = 1
    else:
        iq1 = 1
        iq2 = -1
    return pese1, pese2, iq1, iq2


@devicejit
def set_pair_angle(eig):
    if iprdst > 0:
        if iprdst == 4:
            rtest = randomset()
            # gbeta = (1-rmt2/eig)**8
            gbeta = pese1/(pese1+10)
            iprdst_use = 1 if rtest < gbeta else 4
        elif iprdst == 2 and eig < BHPAIR )
            iprdst_use = 1
        else:
            iprdst_use = iprdst

        for ichrg in (1,2):
            if ichrg == 1:
                ese=pese1
            else:
                ese=ese2
                if iprdst == 4:
                    gbeta = ese/(ese+10)
                    rtest = randomset()
                    if rtest < gbeta:
                        iprdst_use = 1
                    else:
                        iprdst_use = 4
            if iprdst_use == 1:
                pse=sqrt(max(0.0,(ese-REST_MASS)*(ese+REST_MASS)))
                COSTHE = randomset()
                COSTHE=1.0-2.0*COSTHE
                SINTHE=REST_MASS*SQRT((1.0-COSTHE)*(1.0+COSTHE))/(PSE*COSTHE+ESE)
                COSTHE=(ESE*COSTHE+PSE)/(PSE*COSTHE+ESE)
            elif iprdst_use == 2:
                # ZBRANG=( (1/111)*Zeff**(1/3) )**2
                ZTARG=ZBRANG[medium]
                # TTEIG=TOTAL INITIAL PHOTON ENERGY IN ELECTRON REST MASS UNITS
                TTEIG=EIG/REST_MASS
                # TTESE=TOTAL FINAL ELECTRON ENERGY IN ELECTRON REST MASS UNITS
                TTESE=ESE/REST_MASS
                # TTPSE=TOTAL FINAL ELECTRON MOMENTUM IN REST_MASS UNITS
                TTPSE=sqrt((TTESE-1.0)*(TTESE+1.0))
                # THIS IS THE RATIO (r IN PIRS0287)
                ESEDEI=TTESE/(TTEIG-TTESE)
                ESEDER=1.0/ESEDEI
                # DETERMINE THE NORMALIZATION
                XIMIN=1.0/(1.0+(3.141593*TTESE)**2)
                # --- Inline replace: $ SET_PAIR_REJECTION_FUNCTION(REJMIN,XIMIN); -----
                if set_pair_rejection_function:
                    <XXX> = set_pair_rejection_function(REJMIN, XIMIN)
                else:
                    REJMIN = 2.0+3.0*(ESEDEI+ESEDER) -
                            4.00*(ESEDEI+ESEDER+1.0-4.0*(XIMIN-0.5)**2)*(
                                1.0+0.25*LOG(
                                    ((1.0+ESEDER)*(1.0+ESEDEI)/(2.*TTEIG))**2+ZTARG*XIMIN**2
                                    )
                                )
                # End inline replace: $ SET_PAIR_REJECTION_FUNCTION(REJMIN,XIMIN); ----
                YA=(2.0/TTEIG)**2
                XITRY=max(0.01,MAX(XIMIN,min(0.5,SQRT(YA/ZTARG))))
                GALPHA=1.0+0.25*log(YA+ZTARG*XITRY**2)
                GBETA=0.5*ZTARG*XITRY/(YA+ZTARG*XITRY**2)
                GALPHA=GALPHA-GBETA*(XITRY-0.5)
                XIMID=GALPHA/(3.0*GBETA)
                if GALPHA >= 0.0:
                    XIMID=0.5-XIMID+SQRT(XIMID**2+0.25)
                else:
                    XIMID=0.5-XIMID-SQRT(XIMID**2+0.25)

                XIMID=max(0.01,MAX(XIMIN,min(0.5,XIMID)))
                # --- Inline replace: $ SET_PAIR_REJECTION_FUNCTION(REJMID,XIMID); -----
                if set_pair_rejection_function:
                    <XXX> = set_pair_rejection_function(REJMID, XIMID)
                else:
                    REJMID = 2.0+3.0*(ESEDEI+ESEDER) -
                            4.00*(ESEDEI+ESEDER+1.0-4.0*(XIMID-0.5)**2)*(
                                1.0+0.25*LOG(
                                    ((1.0+ESEDER)*(1.0+ESEDEI)/(2.*TTEIG))**2+ZTARG*XIMID**2
                                    )
                                )

                # End inline replace: $ SET_PAIR_REJECTION_FUNCTION(REJMID,XIMID); ----
                # ESTIMATE MAXIMUM OF THE REJECTION FUNCTION
                # FOR LATER USE BY THE REJECTION TECHNIQUE
                REJTOP=1.02*max(REJMIN,REJMID)
                while True:
                    XITST = randomset()
                    # --- Inline replace: $ SET_PAIR_REJECTION_FUNCTION(REJTST,XITST); -----
                    if set_pair_rejection_function:
                        <XXX> = set_pair_rejection_function(REJTST, XITST)
                    else:
                        REJTST = 2.0+3.0*(ESEDEI+ESEDER) -
                                4.00*(ESEDEI+ESEDER+1.0-4.0*(XITST-0.5)**2)*(
                                    1.0+0.25*LOG(
                                        ((1.0+ESEDER)*(1.0+ESEDEI)/(2.*TTEIG))**2+ZTARG*XITST**2
                                        )
                                    )

                    # End inline replace: $ SET_PAIR_REJECTION_FUNCTION(REJTST,XITST); ----
                    RTEST = randomset()
                    # CONVERT THE SUCCESSFUL CANDIDATE XITST TO AN ANGLE
                    THETA=SQRT(1.0/XITST-1.0)/TTESE
                    # LOOP UNTIL REJECTION TECHNIQUE ACCEPTS XITST
                    REJTST_on_REJTOP   = REJTST/REJTOP
                    if RTEST <= REJTST_on_REJTOP and THETA < PI:
                        break
                SINTHE=SIN(THETA);COSTHE=COS(THETA)
            elif iprdst_use == 3:
                COSTHE = randomset()
                COSTHE=1.0-2.0*COSTHE
                sinthe=(1-costhe)*(1+costhe)
                if sinthe > 0:
                    sinthe = sqrt(sinthe)
                else:
                    sinthe = 0
            else:
                # PSE=SQRT(max(1e-10,(ESE-REST_MASS)*(ESE+REST_MASS)))
                # $ RANDOMSET costhe
                # costhe=(ese-(ese+pse)*exp(-2*costhe*log((ese+pse)/REST_MASS)))/pse
                costhe = randomset()
                costhe=1-2*sqrt(costhe)
                sinthe=(1-costhe)*(1+costhe)
                if sinthe > 0:
                    sinthe=sqrt(sinthe)
                else:
                    sinthe=0
            if ichrg == 1:
                UPHI(2,1)
            else:
                sinthe=-sinthe
                NP=NP+1
                UPHI(3,2)
        iq[np] = iq2
        iq(np-1) = iq1
        return
    else:
        THETA=0  # THETA=REST_MASS/EIG
    return iq1, iq2, costhe, sinthe ??, ......??



# ******************************************************************
#                                National Research Council of Canada
@devicejit
def PAIR(rng_states, gid, p):
#
# ******************************************************************
#    For a photon energy below 2.1 MeV, the energies of the pair
#    particles are uniformly distributed in the allowed range via
#    the default replacement for $ SELECT-LOW-ENERGY-PAIR-PRODICTION
#    If the user has a better approach, modify this macro.
#    For a photon energy between 2.1 and 50 MeV the Bethe-Heitler
#    cross section is employed, above 50 MeV the Coulomb-corrected
#    Bethe-Heitler is used.
#    Modified from its original version to make compatible with the
#    changes made in BREMS.
#
#    I. Kawrakow
# ******************************************************************


    NPold = NP  # Set the old stack counter
    medium = p.region.medium
    if i_play_RR == 1:
        #  The user wants to play Russian Roulette. For pair
        #  it is much more efficient to do it BEFORE the
        #  actual sampling
        i_survived_RR = 0  # flag they all survive inititally
        if prob_RR <= 0:
            if n_RR_warning < MAX_RR_WARNING:
                n_RR_warning = n_RR_warning + 1
                logger.warning("Attempt to play Russian Roulette with prob_RR<0! ")
        else:
            rnno_RR = randomset()
            if rnno_RR > prob_RR:  # The pair was killed
                i_survived_RR = 2  # flag both particles eliminated
                if np > 1:
                    np = np-1
                else:
                    #  get a proper exit from PHOTO, we have to leave at least
                    #  one particle on the stack
                    wt[np] = 0
                    e[np] = 0
                return
            else:
                wt[np] = wt[np] / prob_RR

    PEIG=e[np]  # PRECISE ENERGY OF INCIDENT GAMMA
    EIG=PEIG  # ENERGY OF INCIDENT GAMMA
    do_nrc_pair = False

    if itriplet > 0 and eig > 4*REST_MASS:
        itrip = dli_triplet*gle + bli_triplet
        ftrip = medium.a_triplet[itrip]*gle + medium.b_triplet[itrip]
        rnno34 = randomset()
        if rnno34 < ftrip:
            #  Triplet production
            sample_triplet(rng_states, gid, p)
            return


    if pair_nrc == 1:
        # Sample from the NRC pair cross section data base
        # (privided the energy is within the available range)
        k = eig/REST_MASS
        if k < nrcp_emax:
            do_nrc_pair = True
            if k <= nrcp_emin:
                ibin = 1;
            else:
                abin = 1 + log((k-2)/(nrcp_emin-2))*nrcp_dlei
                ibin = abin
                abin = abin - ibin
                rbin = randomset()
                if rbin < abin:
                    ibin = ibin + 1

            xx = alias_sample1(NRC_PAIR_NX_1,nrcp_xdata,
                    nrcp_fdata(1,ibin,medium),nrcp_wdata(1,ibin,medium),
                    nrcp_idata(1,ibin,medium))
            #  The above returns the energy fraction of the positron
            if xx > 0.5:
                pese1 = prm*(1 + xx*(k-2))
                iq1 = 1
                pese2 = peig - pese1
                iq2 = -1
            else:
                pese2 = prm*(1 + xx*(k-2))
                iq2 = 1
                pese1 = peig - pese2
                iq1 = -1

    if not do_nrc_pair:
        if EIG <= 2.1:
            #    BELOW 2.1,USE APPROXIMATION
            # --- Inline replace: $ SELECT_LOW_ENERGY_PAIR_PRODICTION; -----
            pese1, pese2, iq1, iq2 = select_low_energy_pair_production(peig)
            # End inline replace: $ SELECT_LOW_ENERGY_PAIR_PRODICTION; ----
    else:
        #    DECIDE WHETHER TO USE BETHE-HEITLER or BH
        #    COULOMB CORRECTED
        if EIG < 50.0:
            # Use BH without Coulomb correction
            L = 5
            L1 = L + 1

            # Find the actual rejection maximum for this photon energy
            delta = 4*delcm[medium]/eig
            if delta < 1:
                Amax = dl1(l,medium)+delta*(dl2(l,medium)+delta*dl3(l,medium))
                Bmax = dl1(l1,medium)+delta*(dl2(l1,medium)+delta*dl3(l1,medium))
            else:
                aux2 = log(delta+dl6(l,medium))
                Amax = dl4(l,medium)+dl5(l,medium)*aux2
                Bmax = dl4(l1,medium)+dl5(l1,medium)*aux2

            # and then calculate the probability for sampling from (br-1/2)**2
            aux1 = 1 - rmt2/eig
            aux1 = aux1*aux1
            aux1 = aux1*Amax/3
            aux1 = aux1/(Bmax+aux1)
        else:
            # Use BH Coulomb-corrected
            L = 7
            # The absolute maxima are close to the actual maxima at high energies
            # =>use the absolute maxima to save time
            Amax = dl1(l,medium)
            Bmax = dl1(l+1,medium)
            aux1 = bpar(2,medium)*(1-bpar(1,medium)*REST_MASS/eig)

        del0 = eig*delcm[medium]
        Eavail = eig - rmt2

        while True:
            rnno30 = randomset()
            rnno31 = randomset()
            rnno34 = randomset()
            if rnno30 > aux1:  # use the uniform part
                br = 0.5*rnno31
                rejmax = Bmax
                l1 = l+1
            else:  # use the (br-1/2)**2 part of the distribution
                rnno32 = randomset() rnno33 = randomset()
                br = 0.5*(1-max(rnno31,rnno32,rnno33))
                rejmax = Amax
                l1 = l

            Eminus = br*Eavail + REST_MASS
            Eplus  = eig - Eminus
            delta = del0/(Eminus*Eplus)
            if delta < 1:
                rejf = dl1(l1,medium)+delta*(dl2(l1,medium)+delta*dl3(l1,medium))
            else:
                rejf = dl4(l1,medium)+dl5(l1,medium)*log(delta+dl6(l1,medium))

            if rnno34*rejmax <= rejf:
                break

        pese2 = Eminus
        pese1 = peig - pese2
        rnno34 = randomset()
        if rnno34 < 0.5:
            iq1 = -1
            iq2 = 1
        else:
            iq1 = 1
            iq2 = -1


    #    ENERGY GOING TO LOWER SECONDARY HAS NOW BEEN DETERMINED
    ESE2=PESE2
    e[np]=PESE1
    E(NP+1)=PESE2
    #    THIS AVERAGE ANGLE OF EMISSION FOR BOTH PAIR PRODUCTION AND
    #    BREMSSTRAHLUNG IS MUCH SMALLER THAN THE AVERAGE ANGLE OF
    #    MULTIPLE SCATTERING FOR DELTA T TRANSPORT=0.01 R.L.
    #    THE INITIAL AND FINAL MOMENTA ARE COPLANAR
    #    SET UP A NEW 'ELECTRON'
    # --- Inline replace: $ SET_PAIR_ANGLE; -----
    iq1, iq2, costhe, sinthe ??, ......?? = set_pair_angle(eig, ...?)

    # End inline replace: $ SET_PAIR_ANGLE; ----

    #  DEFAULT FOR $ SET-PAIR-ANGLE; is to select the angle from the leading term
    #  of the angular distribution
    CALL UPHI(1,1)
    #    SET UP A NEW 'ELECTRON'
    NP=NP+1
    SINTHE=-SINTHE
    CALL UPHI(3,2)
    iq[np]=iq2
    IQ(NP-1)=iq1
