# angles.py
from typing import Tuple
from math import sqrt

from egsnrc.commons import *
from egsnrc.params import *
from egsnrc.randoms import randomset

import logging
logger = logging.getLogger("egsnrc")

def select_azimuthal_angle() -> Tuple[float, float]:
    """Azimuthal angle selection using a sampling within a box method

    From egsnrc.macros, $SELECT-AZIMUTHAL-ANGLE
    Choose a point randomly within a box such that
       -1 <= x <= 1 and 0 <= y < = 1
    Reject the set if it lies without the inscribed unit semicircle centered
       at (x,y) = (0,0)
    once out of the loop, use the trigonimetric relations (TeX notation)
        \cos 2\phi = (x^2 - y^2)/(x^2 + y^2)
        \sin 2\phi = 2xy/(x^2 + y^2)

    Returns
    -------
    float, float
        cosphi, sinphi
    """

    rhophi2: float = 2  # dummy to start the loop
    while rhophi2 > 1:  # LOOP
        xphi:float = randomset()
        # logger.debug(f"xphi random: {xphi}")
        xphi  = 2*xphi - 1
        xphi2: float = xphi*xphi
        yphi: float = randomset()
        yphi2: float  = yphi * yphi
        rhophi2 = xphi2 + yphi2

    rhophi2 = 1 / rhophi2
    return (xphi2 - yphi2)*rhophi2, 2*xphi*yphi*rhophi2


# Save for connected calls to uphi. See comments below
a: float = 0.0
b: float = 0.0
c: float = 0.0


def uphi(ientry: int, lvl: int):
    """Uniform phi distribution

    Set coordinates for new particle or reset direction cosines of
    old one. Generate random azimuth selection and replace the
    direction cosines with their new values.

    Parameters
    ----------
    ientry: int
        Entry point for the function
        1 for uphi
        2 for uphi2
        3 for the second of two particles when we know two particles
          have a relationship in their corrections.

    lvl: int
        lvl=1 -- old particle, save its direction and adjust it
        lvl=2 -- new particle. adjust direction using saved a,b,c
        lvl=3 -- bremsstrahlung gamma.  save electron direction (next
                 to top of stack), and then adjust gamma direction.

    Note
    ----
        `sinthe` and `costhe` can be changed outside through common.
        `lvl` is a parameter telling which particles to work with.
        Theta (`sinthe` and `costhe`) are always relative to the direction
        of the incident particle before its direction was adjusted.
        Thus when two particles need to have their directions computed,
        the original incident direction is saved in the variable `a`,`b`,`c`
        so that it can be used on both calls."
    """
    global a, b, c

    # "local variables"
    # $real
    # cthet:  5/2*pi-theta, used to evaluate cos(theta) using the sine table
    # rnno38: random number for azimuthal angle selection
    # phi:    azimuthal scattering angle
    # cphi:   5/2*pi-phi
    # a,b,c:  direction cosines before rotation
    # sinps2: sinps2 = a*a + b*b
    # sinpsi: sqrt(sinps2)
    # us,vs:  x- and y- component of scattering vector
    # sindel,cosdel: aux. variables for the rotation algorithm

    # $integer
    # iarg:   index for ausgab
    # lphi,ltheta,lcthet,lcphi:  indices for sine table

    # $define-variables-for-select-azimuthal-angle;
    # save cthet,phi,cphi,a,b,c,sinps2,sinpsi,us,vs,sindel,cosdel;

    np_m1 = np - 1  # ** 0-based
    if iausfl[UPHIAUSB-1+1] != 0:  # ** 0-based
        ausgab(UPHIAUSB)

    assert ientry in (1, 2, 3), f"Invalid UPHI ientry={ientry}"

    # go to (:uphi:,:uphi2:,:nrk:),ientry;
    if ientry == 1:
        raise NotImplementedError(f"Have not coded UPHI for ientry={ientry}")

    # :uphi:; "note: afb 88/12/12 added semi-colon, else bug when overriding sin"
    #         "table look-up"
    # $set interval theta,sinc;
    # $evaluate sinthe using sin(theta);
    # cthet=pi5d2-theta;$set interval cthet,sinc;
    # $evaluate costhe using sin(cthet);

    # Use the following entry if sinthe and costhe are already known.
    # Select phi uniformly over the interval (0,two pi). Then use
    #    pwlf of sin function to get sin(phi) and cos(phi).  The cosine
    #    is gotten by cos(phi)=sin(9*pi/4 - phi).

    # :uphi2:;
    if ientry == 2:
        # it is much faster to use the box method for azimuthal angle selection"
        # than the following                                                   "
        # $randomset rnno38;
        # phi=rnno38*twopi;$set interval phi,sinc;
        # $evaluate sinphi using sin(phi);
        # cphi=pi5d2-phi;$set interval cphi,sinc;
        # $evaluate cosphi using sin(cphi);

        # $select-azimuthal-angle(cosphi,sinphi);

        uphiot.cosphi, uphiot.sinphi = select_azimuthal_angle()

    # :nrk:
    # go to (:old-particle:,:new-particle:,:brems-gamma:),lvl;

    if lvl == 1:  # old particle
        a, b, c = u[np_m1], v[np_m1], w[np_m1]
        # go to adjust
    else:  # lvl 2 or 3;  2->new particle;  3->brems gamma
        np_m2 = np_m1 - 1
        if lvl == 3:  # brems gamma
            raise NotImplementedError(f"Have not coded UPHI for lvl={lvl}")
            # :brems-gamma:
            # a = u[np_m2]
            # b = v[np_m2]
            # c = w[np-m2]
        elif lvl != 2:
            raise ValueError(
                f' STOPPED IN UPHI WITH IENTRY,LVL={ientry},{lvl}'
            )
        # transfer properties to (np) FROM (np-1)
        x[np_m1] = x[np_m2]
        y[np_m1] = y[np_m2]
        z[np_m1] = z[np_m2]
        ir[np_m1] = ir[np_m2]
        wt[np_m1] = wt[np_m2]
        dnear[np_m1] = dnear[np_m2]
        latch[np_m1] = latch[np_m2]

    #   See H.H. Nagel dissertation for coordinate system description.
    #   A rotation is performed to transform direction cosines of the
    #   particle back to the physical frame (from the transport frame)

    # adjust:
    sinps2 = a*a + b*b
    # if sinps2 is small, no rotation is needed
    if sinps2 < 1.0e-20:   # small polar angle case
        u[np_m1] = sinthe * cosphi
        v[np_m1] = sinthe * sinphi
        w[np_m1] = c * costhe    # fixed march 2001 from =costhe"
    else:  # large polar angle case
        sinpsi = sqrt(sinps2)
        us = sinthe * cosphi
        vs = sinthe * sinphi
        sindel = b / sinpsi
        cosdel = a / sinpsi
        u[np_m1] = c*cosdel*us - sindel*vs + a*costhe
        v[np_m1] = c*sindel*us + cosdel*vs + b*costhe
        w[np_m1] = -sinpsi*us + c*costhe

    # $auscall($uphiausa);
    if iausfl[UPHIAUSA-1+1] != 0:  # ** 0-based
        ausgab(UPHIAUSA)
