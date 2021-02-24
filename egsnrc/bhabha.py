from egsnrc.commons import *
from egsnrc.params import *
from egsnrc.randoms import randomset
from .angles import uphi
from math import sqrt

# EMPTY CALLBACKS ----
check_stack = None

#                                National Research Council of Canada
def bhabha():
    """Bhabha scattering of sufficient energy to transport discretely

    Discrete Bhabha scattering (a call to this routine) has been
    arbitrarily defined and calculated to mean Bhabha scatterings
    which impart to the secondary electron sufficient energy that
    it be transported discretely, i.e. e=ae or t=te.  It is not
    guaranteed that the final positron will have this much energy
    however. The exact Bhabha differential cross section is used.
    """

    # $ comin_bhabha # DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
    # COMIN/DEBUG,EGS-VARIANCE-REDUCTION,STACK,
    # THRESH,UPHIOT,USEFUL,RANDOM/
    global np

    stack.npold = np  # Set the old stack counter
    np_m1 = np - 1
    medium_m1 = medium - 1
    peip = e[np_m1]  # PRECISE ENERGY OF INCIDENT POSITRON
    eip = peip  # ENERGY OF INCIDENT POSITRON
    pekin = peip - prm  # PRECISE K.E. OF INCIDENT POSITRON
    ekin = pekin
    t0 = ekin / rm
    e0 = t0 + 1.
    yy = 1. / (t0 + 2.)
    e02 = e0 * e0
    # BETAI2 = E02 / (E02 - 1.)  BLIF 96/2/1 -- not needed for Bhabha fix-up
    beta2 = (e02 - 1.) / e02  # BLIF 96/2/1 -- needed for Bhabha fix-up
    ep0 = te[medium_m1] / ekin
    ep0c = 1. - ep0
    y2 = yy * yy
    yp = 1. - 2. * yy
    yp2 = yp * yp
    b4 = yp2 * yp
    b3 = b4 + yp2
    b2 = yp * (3. + y2)
    b1 = 2. - y2
    #    SAMPLE BR FROM MINIMUM(EP0) TO 1.
    while True:
        br = ep0 / (1. - ep0c * randomset())
        # Apply rejection function.  BLIF 96/2/1 -- Bhabha fix-up
        # rejf2 = ep0c * (betai2-br*(b1-br*(b2-br*(b3-br*b4))))
        rejf2 = (1.0 - beta2 * br * (b1 - br * (b2 - br * (b3 - br * b4))))
        if randomset() <= rejf2:
            break  # leave loop

    # If e- got more than e+, move the e+ pointer and reflect b
    if np + 1 > MXSTACK:
        raise OverflowError(
            f'\n\n In subroutine BHABHA stack size exceeded!\n'
            f' $MAXSTACK = {MXSTACK:9d}, np = {np + 1:9d}'
        )

    if br < 0.5:
        iq[np_m1+1] = -1
    else:
        iq[np_m1] = -1
        iq[np_m1+1] = 1
        br = 1. - br
    # The above puts e+ on top of stack if it has less energy
    #    divide up the energy
    br = max(br, 0.0)  # avoids possible negative number due to round-off
    pekse2 = br * ekin  # precise kinetic energy of secondary 'electron' 2
    pese1 = peip - pekse2  # precise energy of secondary 'electron' 1
    pese2 = pekse2 + prm  # precise energy of secondary 'electron' 2
    ese1 = pese1
    ese2 = pese2
    e[np_m1] = pese1
    e[np_m1+1] = pese2
    #    Bhabha angles are uniquely determined by kinematics
    h1 = (peip + prm) / pekin
    #    direction cosine change for 'old' electron

    # AFB modified the following statement 92/10/28 to avoid
    # numerical difficulties
    # dcosth = h1 * (pese1 - prm) / (pese1 + prm)
    dcosth = min(1.0, h1 * (pese1 - prm) / (pese1 + prm))

    uphiot.sinthe = sqrt(1. - dcosth)
    uphiot.costhe = sqrt(dcosth)
    uphi(2, 1)
    stack.np += 1
    dcosth = h1 * (pese2 - prm) / (pese2 + prm)
    uphiot.sinthe = -sqrt(1. - dcosth)
    uphiot.costhe = sqrt(dcosth)
    uphi(3, 2)