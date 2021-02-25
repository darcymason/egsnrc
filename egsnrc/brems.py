from math import log, sqrt, exp
from .angles import select_azimuthal_angle
from .commons import *
from .params import *
from egsnrc.randoms import randomset

#                                National Research Council of Canada
def brems():
    """Samples bremsstrahlung energy

    Samples bremsstrahlung energy using
     - Coulomb corrected Bethe-Heitler above 50 MeV
     - Bethe-Heitler below 50 MeV
    if ibr_nist == 0, or
     - the NIST bremsstrahlung cross section data base
       (prepared in a form of an alias table for rapid sampling)
    if ibr_nist == 1  or
     - the NRC bremsstrahlung cross section data base, which is
       the same as the NIST database, but with corrections to
       the electron-electron contribution, which are mostly
       important for low Z and low k
    if ibr_nist == 2
    and direction using
     - formula 2BS from from Koch and Motz if IBRDST=1
     - leading term of the brems angular dsstr. if IBRDST=0
     - photon direction = electron direction if IBRDST<0

    This version replaces the original EGS4 implementation
    because of a bug discovered in the EGS4 brems routine
    In order to work properly, the parameter DL1,..,DL6
    are re-calculated in subroutine fix_brems which is called
    from HATCH
    In addition, this version has the internal capability of
    bremsstrahlung splitting.
    To use bremsstrahlung splitting, set nbr_split (COMON/BREMPR/)
    to the desired number > 1 (1 is the default)
    Be aware that event-by-event energy conservation is NOT
    guaranteed, so don't use for calculations where this is
    important (e.g. calculation of detector response functions)
    The result will be nbr_split photons, all with the weight
    wt(npold)/nbr_split, and an electron with the original weight
    and energy given by the incident energy - energy of last photon

    I. Kawrakow, January 2000

    D. Mason, Feb 2021 - transpile to Python
    """
    # $ comin_brems # DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
    # COMIN/DEBUG,BREMPR,EGS-VARIANCE-REDUCTION,
    # STACK,THRESH,UPHIOT,USEFUL,RANDOM/
    # and common/nist_brems/  - for log_ap  DLM 2021-02

    # $ define_local_variables_brems
    # real*8 z2max,z2maxi,aux1,aux3,aux4,aux5,aux2,weight

    # $ENERGY PRECISION
    #   peie,   precise incident electron energy
    #   pesg,   presice energy of emitted photon
    #   pese;   precise total energy of scattered electron
    # $REAL
    #   eie,    total incident electron energy
    #   ekin,   kinetic incident energy
    #   brmin,  ap[medium] / ekin
    #   waux,   for faster sampling of 1/br
    #   aux,    ese / eie
    #   ajj,    for energy bin determination if alias sampling is employed
    #   alias_sample1,
    #   br,     energy fraction of secondary photon
    #   esg,    energy of secondary photon
    #   ese,    total energy of secondary electron
    #   delta,  scaled momentum transfer
    #   phi1,   screening function
    #   phi2,   screening function
    #   rejf;   screening rejection function
    #
    # Brems angle selection variables
    # $REAL
    #   a,b,c,  direction cosines of incident `electron'
    #   sinpsi, sindel, cosdel, us, vs,
    #           all used for rotations
    #   ztarg,  (Zeff**1/3/111)**2, used for 2BS angle sampling
    #   tteie,  total energy in units of rest energy
    #   beta,   electron velocity in units of speed of light
    #   y2max,  maximum possible scaled angle
    #   y2maxi, inverse of the above
    #   ttese,  new electron energy in units of rm
    #   rjarg1,rjarg2,rjarg3,rejmin,rejmid,rejmax,rejtop,rejtst,
    #           all of them used for angle rejection function calcs
    #   esedei, new total energy over old total energy
    #   y2tst,  scaled angle, costhe = 1 - 2*y2tst/y2max
    #   y2tst1,
    #   rtest,  random number for rejection
    #   xphi,yphi,xphi2,yphi2,rhophi2,cphi,sphi;
    #           all of the above is for azimuthal angle sampling
    #
    # $INTEGER
    #   l, l1, ibr, jj, j

    # The user can turn off brems production by setting nbr_split to zero
    if nbr_split < 1:
        return

    stack.npold = np  # Set the old stack counter

    # ** 0-based Python
    np_m1 = np - 1
    npold_m1 = npold - 1
    medium_m1 = medium - 1

    peie = e[np_m1]  # Precise energy of incident 'electron'
    eie = peie  # Energy of incident 'electron'
    weight = wt[np_m1] / nbr_split

    # Decide which distribution to use (b-h coulomb corrected is
    # used from 50 to 20000 mev, b-h is used 1.5 to 50 mev)

    l = 1 if eie < 50.0 else 3
    l1 = l + 1
    l_m1 = l - 1  # ** 0-based Python
    l1_m1 = l1 - 1

    ekin = peie - prm
    brmin = ap[medium_m1] / ekin
    # waux = -log(brmin)
    waux = elke - log_ap[medium_m1] # this saves the time consuming log evaluation
                                # log_ap = log(ap[medium_m1]) is calculated in
                                # fix_brems for each medium, elke is needed
                                # in electr to calculate the branching ratios
                                # and therefore it must be known at this point

    if ibrdst >= 0:
        # inrdst >=0 means we will sample the photon emmision angle
        # from KM-2BS (ibrdst=1) or from the leading term (ibrdst=0).
        # If nbr_split > 1, we can re-use the following quantities several time
        a = u[np_m1]
        b = v[np_m1]
        c = w[np_m1]
        sinpsi = a*a + b*b
        if sinpsi > 1e-20:
            sinpsi = sqrt(sinpsi)
            sindel = b / sinpsi
            cosdel = a / sinpsi

        ztarg = zbrang[medium_m1]
        tteie = eie / rm
        beta = sqrt((tteie - 1)*(tteie + 1)) / tteie
        y2max = 2*beta*(1 + beta)*tteie*tteie
        y2maxi = 1 / y2max
        if ibrdst == 1:
            z2max = y2max + 1
            z2maxi = sqrt(z2max)

    if ibr_nist >= 1:
        ajj = 1 + (
            waux + log_ap[medium_m1] - nb_lemin[medium_m1]
            ) * nb_dlei[medium_m1]
        jj = ajj
        ajj = ajj - jj
        if jj > MXBRES:
            jj = MXBRES
            ajj = -1

    for ibr in range(nbr_split):
        if ibr_nist >= 1:
            # use the NIST or NRC bremsstrahlung cross section data base
            if ekin > nb_emin[medium_m1]:
                j = jj + 1 if randomset() < ajj else jj
                j_m1 = j - 1  # ** 0-based
                br = egsfortran.alias_sample1(
                    MXBRXS,
                    nb_xdata[0, j_m1, medium_m1],
                    nb_fdata[0, j_m1, medium_m1],
                    nb_wdata[1, j_m1, medium_m1],
                    nb_idata[1, j_m1, medium_m1]
                )
            else:
                br = randomset()
            esg = ap[medium_m1]*exp(br*waux)
            pesg = esg
            pese = peie - pesg
            ese = pese
        else:
            while True:  # User wants to use Bethe-Heitler
                br = brmin * exp(randomset() * waux)
                esg = ekin*br
                pesg = esg
                pese = peie - pesg
                ese = pese
                delta = esg / eie / ese * delcm[medium_m1]
                aux = ese / eie
                if delta < 1:
                    phi1 = (
                        dl1[l_m1, medium_m1]
                        + delta * (
                            dl2[l_m1, medium_m1] + delta*dl3[l_m1, medium_m1]
                        )
                    )
                    phi2 = (
                        dl1[l1_m1, medium_m1]
                        + delta * (
                            dl2[l1_m1, medium_m1] + delta*dl3[l1_m1, medium_m1]
                        )
                    )
                else:
                    phi1 = (
                        dl4[l_m1, medium_m1]
                        + dl5[l_m1, medium_m1] * log(
                            delta + dl6[l_m1, medium_m1]
                        )
                    )
                    phi2 = phi1

                rejf = (1 + aux*aux)*phi1 - 2*aux*phi2 / 3
                if randomset() < rejf:
                    break  # exit loop

        # Set up the new photon
        stack.np += 1
        np_m1 = np - 1  # ** 0-based

        if np > MXSTACK:
            raise OverflowError(
                f'\n\n Stack overflow in BREMS! np = {np+1}.'
                ' Increase MXSTACK and try again'
            )

        e[np_m1] = pesg
        iq[np_m1] = 0

        # transfer properties to (np) FROM (np-1)
        np_m2 = np_m1 - 1
        x[np_m1] = x[np_m2]
        y[np_m1] = y[np_m2]
        z[np_m1] = z[np_m2]
        ir[np_m1] = ir[np_m2]
        wt[np_m1] = wt[np_m2]
        dnear[np_m1] = dnear[np_m2]
        latch[np_m1] = latch[np_m2]

        wt[np_m1] = weight
        if ibrdst < 0:
            # The photon will inherit the direction from the electron.
            # This option is given so that the user can implement their own
            # brems angle schemes via a call to ausgab
            u[np_m1] = u[npold_m1]
            v[np_m1] = v[npold_m1]
            w[np_m1] = w[npold_m1]
        else:
            if ibrdst == 1:
                # See after function for original implementation comments
                ttese = ese / rm
                esedei = ttese / tteie
                rjarg1 = 1 + esedei*esedei
                rjarg2 = rjarg1 + 2*esedei
                aux = 2*ese*tteie / esg
                aux = aux*aux
                aux1 = aux*ztarg
                if aux1 > 10:
                    rjarg3 = lzbrang[medium_m1] + (1 - aux1) / aux1**2
                else:
                    rjarg3 = log(aux / (1 + aux1))

                rejmax = rjarg1*rjarg3 - rjarg2
                while True:
                    y2tst = randomset()
                    rtest = randomset()
                    aux3 = z2maxi / (y2tst + (1 - y2tst) * z2maxi)
                    rtest = rtest * aux3 *rejmax
                    y2tst = aux3**2 - 1
                    y2tst1 = esedei * y2tst / aux3**4
                    aux4 = 16 * y2tst1 - rjarg2
                    aux5 = rjarg1 - 4 * y2tst1
                    if rtest < aux4 + aux5 * rjarg3:
                        break
                    aux2 = log(aux / (1 + aux1 / aux3**4))
                    rejtst = aux4 + aux5*aux2
                    if rtest < rejtst:
                        break
            else:
                y2tst = randomset() / (1 - y2tst + y2maxi)

            uphiot.costhe = 1 - 2 * y2tst * y2maxi
            uphiot.sinthe = sqrt(max((1 - costhe)*(1 + costhe), 0.0))
            # --- Inline replace: $ SELECT_AZIMUTHAL_ANGLE(cphi,sphi); -----
            cphi, sphi = select_azimuthal_angle()
            if sinpsi >= 1e-10:
                us = sinthe*cphi
                vs = sinthe*sphi
                u[np_m1] = c*cosdel*us - sindel*vs + a*costhe
                v[np_m1] = c*sindel*us + cosdel*vs + b*costhe
                w[np_m1] = c*costhe - sinpsi*us
            else:
                u[np_m1] = sinthe*cphi
                v[np_m1] = sinthe*sphi
                w[np_m1] = c*costhe

    e[npold_m1] = pese


#             This is the original implementation [under `if ibrdst == 1:`]
#             suggested by Alex Bielajew. Commented out as
#             the implementation above is way more efficient.
#             IK, Sep. 2004.
#         ttese = ese / rm
#         esedei = ttese / tteie
#         rjarg1 = 1 + esedei*esedei
#         rjarg2 = 3*rjarg1 - 2*esedei
#         rjarg3 = ((1 - esedei) / (2*tteie*esedei))**2

# ; Y2TST1 = (1. + 0.0)**2
# REJMIN= (4. + log(RJARG3 + ZTARG / Y2TST1))*(4.*ESEDEI*0.0 / Y2TST1 - RJARG1) + RJARG2


# ; Y2TST1 = (1. + 1.0)**2
# REJMID= (4. + log(RJARG3 + ZTARG / Y2TST1))*(4.*ESEDEI*1.0 / Y2TST1 - RJARG1) + RJARG2


# ; Y2TST1 = (1. + y2max)**2
# REJMAX= (4. + log(RJARG3 + ZTARG / Y2TST1))*(4.*ESEDEI*y2max / Y2TST1 - RJARG1) + RJARG2

#         rejtop = max(rejmin,rejmid,rejmax)
#         LOOP [
#             y2tst = randomset() y2tst = y2tst / (1 - y2tst + y2maxi)

# ; Y2TST1 = (1. + Y2TST)**2
# REJTST= (4. + log(RJARG3 + ZTARG / Y2TST1))*(4.*ESEDEI*Y2TST / Y2TST1 - RJARG1) + RJARG2

#             rtest = randomset()
#         if rtest*rejtop <= REJTST:
#             break
#         */