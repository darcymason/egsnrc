from math import sqrt
from .commons import *
from .params import *
from .angles import uphi
from egsnrc.randoms import randomset


#                                National Research Council of Canada
def moller():
    """Moller scattering of sufficient energy to transport discretely

        Discrete moller scattering (a call to this routine) has been
        arbitrarily defined and calculated to mean moller scatterings
        which impart to the secondary electron sufficient energy that
        it be transported discretely.  The threshold to transport an
        electron discretely is a total energy of ae or a kinetic energy
        of te=ae-rm.  Since the kinetic energy transfer is always, by
        definition, less than half of the incident kinetic energy, this
        implies that the incident energy, eie, must be larger than
        thmoll=te*2 + rm.  The rest of the collision contribution is
        subtracted continuously from the electron as ionization
        loss during transport.
    """
    # COMIN/DEBUG,STACK,THRESH,UPHIOT,USEFUL,RANDOM,EGS-IO,
    # EGS-VARIANCE-REDUCTION/
    # EII-DATA (ei_xxx...), EDGE (binding_energies)  - DLM, 2021-02

    # Local variables
    # $ENERGY PRECISION
    #      peie,   precise total energy of incident electron
    #      pekse2, precise kinetic energy of 2nd secondary electron
    #      pese1,  precise total energy of 1st secondary electron
    #      pese2,  precise total energy of 2nd secondary electron
    #      pekin,  precise kinetic energy of incident electron
    #      h1,     used for polar scattering angle calculation
    #      dcosth; polar scattering angle squared
    # $REAL
    #      eie,    total energy of incident electron
    #      ekin,   kinetic energy of incident electron
    #      t0,     kinetic energy of incident electron in units of rm
    #      e0,     total energy of incident electron in units of RM
    #      extrae, energy above the Moller threshold
    #      e02,    e0**2
    #      ep0,    minimum alowed kinetic energy fraction
    #      g2,g3,  used for rejection function calculation
    #      gmax,   maximum value of the rejection function
    #      br,     kinetic energy fraction to lowew energy electron
    #      r,      (1-br) / br
    #      rejf4,  rejection function
    #      rnno28, random number for rejection
    #      ese1,   energy of 1st secondary electron
    #      ese2;   energy of 2nd secondary electron

    # real*8 sigm, pbrem, rsh, Uj, sig_j
    # integer*4 lelke, iele, ish, nsh, ifirst, i, jj, iZ, iarg

    # IRCODE = 1;  appears to be unused, IK Oct 97
    stack.npold = np # Set the old stack counter

    np_m1 = np - 1  # ** 0-based
    medium_m1 = medium - 1

    peie = e[np_m1] # precise energy of incident electron
    eie = peie # energy of incident electron
    pekin = peie - prm # precise k.e. of incident electron
    ekin = pekin

    if eii_flag > 0 and eii_nsh[medium_m1] > 0:
        # The EII flag is set and this medium has shells for which we want to
        # simulate EII => sample if the interaction is with a EII shell

        # $ SET INTERVAL elke,eke;
        lelke = eke1[medium_m1]*elke + eke0[medium_m1]
        lelke_m1 = lelke - 1  # ** 0-based

        # EVALUATE sigm USING esig(elke)
        sigm = esig1[lelke_m1, medium_m1]*elke+ esig0[lelke_m1, medium_m1]

        # EVALUATE pbrem USING ebr1(elke)
        pbrem = ebr11[lelke_m1, medium_m1]*elke+ ebr10[lelke_m1, medium_m1]

        sigm = sigm * (1 - pbrem)
        rsh = sigm * randomset()
        for iele_m1 in range(nne[medium_m1]):  # ** 0-based iele only as index
            iZ = int(zelem[medium_m1, iele_m1] + 0.5)
            iZ_m1 = iZ -1
            nsh = eii_no[medium_m1, iele_m1]
            if nsh > 0:
                ifirst = eii_first[medium_m1, iele_m1]
                for ish_m1 in range(nsh):
                    Uj = binding_energies[ish_m1, iZ_m1]
                    if ekin > Uj and (Uj > te[medium_m1] or Uj > ap[medium_m1]):
                        jj = ifirst + ish_m1  # orig was ifirst + ish - 1
                        i = eii_a(jj)*elke + eii_b(jj) + (jj - 1)*N_EII_BINS
                        sig_j = eii_xsection_a(i)*elke + eii_xsection_b(i)
                        sig_j = sig_j*pz[medium_m1, iele_m1]*eii_cons[medium_m1]
                        rsh = rsh - sig_j
                        if rsh < 0:
                            if IAUSFL[EIIB -1 + 1] != 0:
                                ausgab(EIIB)
                            egsfortran.eii_sample(ish_m1+1, iZ, Uj)
                            if IAUSFL[EIIA -1 + 1] != 0:
                                ausgab(EIIA)
                            return

    if ekin <= 2*te[medium_m1]:
        return
    t0 = ekin / rm
    e0 = t0 + 1.0
    extrae = eie - thmoll[medium_m1]
    e02 = e0*e0
    # betai2 = e02 / (e02 - 1.0) # BLIF 96/2/1 -- not needed for moller fix-up
    te[medium_m1] / ekin
    # g1 = (1. - 2.*ep0)*betai2 # BLIF 96/2/1 -- not needed for moller fix-up
    g2 = t0*t0 / e02
    g3 = (2.*t0 + 1.) / e02
    #    H.H.Nagel has constructed a factorization of the frequency
    #    distribution function for the moller differential cross
    #    section used as suggested by butcher and messel.
    #    (H.H.Nagel, op.cit., p. 53-55)
    #    However, a much simpler sampling method which does not become
    #    very inefficient near thmoll is the following. . .
    #    Let br = eks / ekin,  where eks is kinetic energy transfered to the
    #    secondary electron and ekin is the incident kinetic energy.

    #    Modified (7 Feb 1974) to use the true moller cross section.
    #    that is, instead of the e+ e- average given in the rossi
    #    formula used by nagel.  The sampling scheme is that
    #    used by Messel and Crawford (EPSDF 1970 p.13)
    #    First sample (1/br**2) over (te/ekin, 1/2) . . .

    gmax = (1. + 1.25*g2)  # BLIF 96/2/1 -- Moller fix-up
    while True:  #  to retry if rejected
        br = te[medium_m1] / (ekin - extrae*randomset())
        #    use Messel and Crawfords rejection function.
        r = br / (1. - br)
        rejf4 = (1. + g2*br*br + r*(r - g3))  # BLIF 96/2/1--remove "g1*..."
        rnno28 = gmax * randomset()  # BLIF 96/2/1 -- Moller fix-up
        if rnno28 < rejf4:   # Try until accepted. end rejection loop
            break  # leave loop

    pekse2 = br*ekin  # precise kinetic energy of secondary electron #2
    pese1 = peie - pekse2  # precise energy of secondary electron #1
    pese2 = pekse2 + prm  # precise energy of secondary electron #2
    e[np_m1] = pese1

    if np + 1 > MXSTACK:
        raise OverflowError(
            f'\n\n In subroutine MOLLER stack size exceeded!\n'
            f' $MAXSTACK = {MXSTACK:9d}, np = {np + 1:9d}'
        )

    e[np_m1 + 1] = pese2
    #    Since `br` < 0.5, `e[np_m1 + 1]` must be < `e[np_m1]`.
    #    Moller angles are uniquely determined by kinematics

    #  One possible way of dealing with double counting of angular
    #  deflections in inelastic scattering would be to
    #  not deflect the 'old' electron as these deflections are
    #  already taken into account in the multiple elastic scattering
    #  This approach has the disadvantage of loosing correlations
    #  between big energy losses and strong angular deflections
    #  The advantage of such an approach is its simplicity.
    #  If spin effects for multiple elastic scattering are turned on,
    #  the double counting is taken into account by the appropriate
    #  modification of the scattering power (which depends on `ae`)
    #
    #
    #  IK, June 1999

    h1 = (peie + prm) / pekin
    # direction cosine change for 'old' electron
    dcosth = h1*(pese1 - prm) / (pese1 + prm)
    uphiot.sinthe = sqrt(1.0 - dcosth)
    uphiot.costhe = sqrt(dcosth)

    # This will turn off the Moller ang. deflections:
    # uphiot.sinthe = 0; uphiot.costhe = 1

    uphi(2,1)

    # Related change and (x,y,z) setup for 'new' electron
    stack.np += 1
    np_m1 = np - 1
    iq[np_m1] = -1
    dcosth = h1*(pese2 - prm) / (pese2 + prm)
    uphiot.sinthe = -sqrt(1.0 - dcosth)
    uphiot.costhe = sqrt(dcosth)
    uphi(3,2)
