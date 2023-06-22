# angles.py
from math import sqrt, sin, cos
from egsnrc.egsrandom import random_kfloat
from egsnrc.constants import TWO_PI
from egsnrc.params import SMALL_POLAR_ANGLE_THRESHOLD
from egsnrc.config import device_jit

@device_jit
def uphi(rng_states, gid, p, sinthe, costhe):
    """Uniform phi distribution

    Returns
    -------
    u, v, w direction cosines
    """
    # unlike original Mortran calculate directly.
    # should be fast enough on todays computers
    phi = TWO_PI * random_kfloat(rng_states, gid)

    sinps2 = p.u * p.u + p.v * p.v
    # Small polar angle case
    if sinps2 < SMALL_POLAR_ANGLE_THRESHOLD:
        u = sinthe * cos(phi)
        v = sinthe * sin(phi)
        w = p.w * costhe
    else:
        sinpsi = sqrt(sinps2)
        us = sinthe * cos(phi)
        vs = sinthe * sin(phi)
        sindel = p.v / sinpsi
        cosdel = p.u / sinpsi

        u = p.w * cosdel * us - sindel * vs + p.u * costhe
        v = p.w * sindel * us + cosdel * vs + p.v * costhe
        w = -sinpsi * us + p.w * costhe
    return u, v, w
