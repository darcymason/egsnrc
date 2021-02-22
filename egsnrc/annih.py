from math import log, sqrt, exp
from .angles import select_azimuthal_angle
from .commons import *
from .params import *
from .randoms import randomset

import logging
logger = logging.getLogger("egsnrc")


def annih():
    """Gamma spectrum for two gamma in-flight positron annihilation.

    Uses scheme based on heitler's p269-27o formulae.

    If the user requests radiative splitting (via nbr_split > 1),
    this routine produces 2*nbr_split annihilation photons at once,
    each carying the fraction 1/nbr_split of the weight of the
    incident positron.

    Except for taking out the calculation of
    log((1.0-EP0)/EP0) out of the sampling loop and using a
    rejection function normalized to its maximum, the sampling
    technique is the same as the original EGS4 implementation.

    I. Kawrakow, January 2000

    Python conversion:  Darcy Mason, Feb 2021
    """

    # $ comin_annih # DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
                    # COMIN/DEBUG,STACK,UPHIOT,USEFUL,RANDOM,
                    # EGS-VARIANCE-REDUCTION/

    # DEFINE-LOCAL-VARIABLES-FOR-ANNIH
    # $ENERGY PRECISION
    #       PAVIP,    "precise total energy in the laboratory frame"
    #       PESG1,    "precise energy of 1st annihilation photon"
    #       PESG2;    "precise energy of 2nd annihilation photon"
    # $REAL AVIP,     "total energy in the laboratory frame"
    #       A,        "total energy in units of the electron's rest energy"
    #       G,T,P,    "energy, kinetic energy and momentum in units of RM"
    #       POT,      "P/T"
    #       EP0,      "minimum fractional energy"
    #       WSAMP,    "to avoid un-necessary calc. of Log((1-ep0)/ep0)"
    #       RNNO01,   "random numbers"
    #       RNNO02,
    #       EP,       "fractional energy of the more energetic photon"
    #       REJF,     "rejection function"
    #       ESG1,     "energy of the more energetic photon"
    #       ESG2,     "energy of the less energetic photon"
    #       aa,bb,cc,sinpsi,sindel,cosdel,us,vs,cphi,sphi;
    #                 "for inline rotations"
    # $INTEGER
    #       ibr;

    stack.npold = np  # Set the old stack counter

    if nbr_split <= 0:
        return

    logger.debug("In annih")
    np_m1 = np - 1  # 0-based
    pavip = e[np_m1] + prm  # PRECISE AVAILABLE ENERGY OF INCIDENT POSITRON,
                            # i.e. electron assumed to be at rest
    avip = pavip  # AVAILABLE ENERGY OF INCIDENT POSITRON
    a = avip / rm
    # AI=1.0/A;  AI not necessary, IK Oct 97
    g = a - 1.0
    t = g - 1.0
    p = sqrt(a*t)
    pot = p / t
    ep0 = 1.0 / (a + p)
    #    SAMPLE 1/EP FROM EP=EP0 TO 1.0-EP0
    # Take the calculation of the logarithm out of the loop, IK Oct 97
    wsamp = log((1.0 - ep0) / ep0)

    aa, bb, cc = u[np_m1], v[np_m1], w[np_m1]
    sinpsi = aa*aa + bb*bb
    if sinpsi > 1e-20:
        sinpsi = sqrt(sinpsi)
        sindel = bb / sinpsi
        cosdel = aa / sinpsi

    if nbr_split > 1:
        wt[np_m1] /= nbr_split

    for ibr in range(nbr_split):  # nbr_split > 1 means we want splitting for any
                                # radiative event
        if np + 1 > MXSTACK:
            raise OverflowError(
                f' Stack overflow in ANNIH! np = {np+1}.'
                ' Increase MXSTACK and try again'
            )

        while True:
            ep = ep0 * exp(randomset() * wsamp)
            #    NOW DECIDE WHETHER TO ACCEPT
            # REJF=1.0-EP+AI*AI*(2.0*G-1.0/EP)
            # The above rejection function has a maximum = 1 - 2/A**2
            # For efficiency, it is better to divide by the maximum value, IK Oct 97
            rejf = 1 - (ep*a - 1)**2 / (ep*(a*a - 2))
            if randomset() <= rejf:
                break

        # SET UP ENERGIES
        esg1 = avip * ep  # ENERGY OF SECONDARY GAMMA 1
        pesg1 = esg1  # PRECISE ENERGY OF SECONDARY GAMMA 1
        e[np_m1] = pesg1
        iq[np_m1] = 0

        ip = npold if ibr == 0 else np - 1  # ibr == 0 due to Python 0-based
        ip_m1 = ip - 1  # ** 0-based arrays
        # transfer properties to (np) FROM (ip)
        x[np_m1] = x[ip_m1]; y[np_m1] = y[ip_m1]; z[np_m1] = z[ip_m1]
        ir[np_m1] = ir[ip_m1]
        wt[np_m1] = wt[ip_m1]
        dnear[np_m1] = dnear[ip_m1]
        latch[np_m1] = latch[ip_m1]

        uphiot.costhe = max(-1.0, min(1.0, (esg1-rm)*pot/esg1))
        uphiot.sinthe = sqrt(1.0 - costhe*costhe)

        cphi, sphi = select_azimuthal_angle()
        if sinpsi >= 1e-10 :
            us = sinthe*cphi
            vs = sinthe*sphi
            u[np_m1] = cc*cosdel*us - sindel*vs + aa*costhe
            v[np_m1] = cc*sindel*us + cosdel*vs + bb*costhe
            w[np_m1] = cc*costhe - sinpsi*us
        else:
            u[np_m1] = sinthe*cphi; v[np_m1] = sinthe*sphi; w[np_m1] = cc*costhe

        stack.np += 1
        np_m1 = np - 1  # 0-based
        pesg2 = pavip - pesg1
        esg2 = pesg2
        e[np_m1] = pesg2
        iq[np_m1] = 0
        np_m2 = np_m1 - 1
        # transfer properties to (np) FROM (np-1)
        x[np_m1] = x[np_m2]; y[np_m1] = y[np_m2]; z[np_m1] = z[np_m2]
        ir[np_m1] = ir[np_m2]
        wt[np_m1] = wt[np_m2]
        dnear[np_m1] = dnear[np_m2]
        latch[np_m1] = latch[np_m2]

        uphiot.costhe = max(-1.0, min(1.0, (esg2-rm)*pot/esg2))
        uphiot.sinthe = -sqrt(1.0 - costhe*costhe)
        if sinpsi >= 1e-10:
            us = sinthe*cphi
            vs = sinthe*sphi
            u[np_m1] = cc*cosdel*us - sindel*vs + aa*costhe
            v[np_m1] = cc*sindel*us + cosdel*vs + bb*costhe
            w[np_m1] = cc*costhe - sinpsi*us
        else:
            u[np_m1] = sinthe*cphi
            v[np_m1] = sinthe*sphi
            w[np_m1] = cc*costhe

        stack.np += 1

    stack.np -= 1

