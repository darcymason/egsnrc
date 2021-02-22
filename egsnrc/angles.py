# angles.py
from typing import Tuple
from egsnrc.randoms import randomset

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
        xphi2:float = xphi*xphi
        yphi:float = randomset()
        yphi2:float  = yphi * yphi
        rhophi2 = xphi2 + yphi2

    rhophi2 = 1 / rhophi2
    return (xphi2 - yphi2)*rhophi2, 2*xphi*yphi*rhophi2
