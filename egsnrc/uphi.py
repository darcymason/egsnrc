from math import sqrt

from egsnrc.commons import *
from egsnrc.params import *
from egsnrc.angles import select_azimuthal_angle

import logging
logger = logging.getLogger("egsnrc")


def uphi(ientry, lvl):
    """uniform phi distribution
    Set coordinates for new particle or reset direction cosines of
    old one. Generate random azimuth selection and replace the
    direction cosines with their new values.
    """
    # "local variables"
    # $real cthet,  "5/2*pi-theta, used to evaluate cos(theta) using the sine table"
    #       rnno38, "random number for azimuthal angle selection"
    #       phi,    "azimuthal scattering angle"
    #       cphi,   "5/2*pi-phi"
    #       a,b,c,  "direction cosines before rotation"
    #       sinps2, "sinps2=a*a+b*b"
    #       sinpsi, "sqrt(sinps2)"
    #       us,vs,  "x- and y- component of scattering vector"
    #       sindel,cosdel;
    #               "aux. variables for the rotation algorithm"

    # $integer
    #       iarg,   "index for ausgab"
    #       lphi,ltheta,lcthet,lcphi;
    #               "indeces for sine table"

    # $define-variables-for-select-azimuthal-angle;
    # save cthet,phi,cphi,a,b,c,sinps2,sinpsi,us,vs,sindel,cosdel;

    np_m1 = np - 1  # ** 0-based
    if iausfl[UPHIAUSB-1+1] != 0:  # ** 0-based
        ausgab(UPHIAUSB)

    # go to (:uphi:,:uphi2:,:nrk:),ientry;
    if ientry != 2:
        raise NotImplementedError(f"Have not coded UPHI for ientry={ientry}")

    # :uphi:; "note: afb 88/12/12 added semi-colon, else bug when overriding sin"
    #         "table look-up"
    # $set interval theta,sinc;
    # $evaluate sinthe using sin(theta);
    # cthet=pi5d2-theta;$set interval cthet,sinc;
    # $evaluate costhe using sin(cthet);

    # "   use the following entry if sinthe and costhe are already known."
    # "   select phi uniformly over the interval (0,two pi). then use    "
    # "   pwlf of sin function to get sin(phi) and cos(phi).  the cosine "
    # "   is gotten by cos(phi)=sin(9*pi/4 - phi).                       "

    # :uphi2:;

    # it is much faster to use the box method for azimuthal angle selection"
    # than the following                                                   "
    # $randomset rnno38;
    # phi=rnno38*twopi;$set interval phi,sinc;
    # $evaluate sinphi using sin(phi);
    # cphi=pi5d2-phi;$set interval cphi,sinc;
    # $evaluate cosphi using sin(cphi);

    # $select-azimuthal-angle(cosphi,sinphi);
    cosphi, sinphi = select_azimuthal_angle()

    #   use the following entry for the second of two particles when we"
    #   know two particles have a relationship in their corrections.   "
    #   note: sinthe and costhe can be changed outside through common. "
    #   lvl is a parameter telling which particles to work with.       "
    #   theta (sinthe and costhe) are always relative to the direction "
    #   of the incident particle before its direction was adjusted.    "
    #   thus when two particles need to have their directions computed,"
    #   the original incident direction is saved in the variable a,b,c "
    #   so that it can be used on both calls."

    #   lvl=1 -- old particle, save its direction and adjust it
    #   lvl=2 -- new particle. adjust direction using saved a,b,c
    #   lvl=3 -- bremsstrahlung gamma.  save electron direction (next
    #   to top of stack), and then adjust gamma direction.

    # :nrk:
    # go to (:old-particle:,:new-particle:,:brems-gamma:),lvl;

    if lvl != 1:
        raise NotImplementedError(f"Have not coded UPHI for lvl={lvl}")

    # :old-particle:

    a, b, c = u[np_m1], v[np_m1], w[np_m1]

    # go to :adjust:;

    # :brems-gamma:
    # a=u(np-1);b=v(np-1);c=w(np-1);

    # :new-particle:
    # $transfer properties to (np) from (np-1);

    #   see h.h. nagel dissertation for coordinate system description.
    #   a rotation is performed to transform direction cosines of the
    #   particle back to the physical frame (from the transport frame)

    # :adjust:
    sinps2 = a*a + b*b
    #   if sinps2 is small, no rotation is needed
    if sinps2 < 1.0e-20:   # small polar angle case
        u[np_m1] = sinthe*cosphi
        v[np_m1] = sinthe*sinphi
        w[np_m1] = c*costhe    # fixed march 2001 from =costhe"
    else:  # large polar angle case
        sinpsi = sqrt(sinps2)
        us = sinthe*cosphi
        vs = sinthe*sinphi
        sindel = b/sinpsi
        cosdel = a/sinpsi
        u[np_m1] = c*cosdel*us-sindel*vs+a*costhe
        v[np_m1] = c*sindel*us+cosdel*vs+b*costhe
        w[np_m1] = -sinpsi*us+c*costhe

    # $auscall($uphiausa);
    if iausfl[UPHIAUSA-1+1] != 0:  # ** 0-based
        ausgab(UPHIAUSA)
