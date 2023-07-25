# angles.py
from math import sqrt, sin, cos
from egsnrc.egsrandom import random_kfloat
from egsnrc.constants import TWO_PI
from egsnrc.params import SMALL_POLAR_ANGLE_THRESHOLD
from egsnrc import config

@config.device_jit
def uphi(rng_states, gid, p, sinthe, costhe):
    """Uniform phi distribution

    Returns
    -------
    u, v, w direction cosines
    """
    # In original Mortran UPHI(ientry, lvl):
    #   lvl only used for ientry==3
    # SAVEd between calls: CTHET,PHI,CPHI,A,B,C,SINPS2,SINPSI,US,VS,SINDEL,COSDEL
    # ientry:
    #   1: calculate sinthe, costhe from theta (falls through to 2)
    #        (only used once in code, in PAIR after SET-PAIR-ANGLE, which
    #         itself calls uphi(2,1) for 1st particle, uphi(3,2) on -sinthe for 2nd),
    #         ?? ** but does not set theta itself??
    #      Fortran code for this step:
    #          SINTHE = sin(THETA)
    #          CTHET  = PI5D2-THETA  # pi5d2 = 5/2*PI, i.e. 2.5*PI
    #          COSTHE = sin(CTHET)
    #   2: Randomly set phi (falls through to 3)
    #   3: Adjust particle directions accordingly:
    #     Lvl
    #         2: new particle: copy n-1 particle, use saved abc, ADJUST
    #
    #         In old code source, only lvl 2 was used, i.e. (3,2).
    #         For reference, others were:
    #         1: old particle: abc = uvw (top particle), then ADJUST
    #         3: brems-gamma: abc=electron(n-1) uvw, ADJUST nth particle
    # :ADJUST:
    #   calc uvw from theta, phi:
    #       "a rotation is performed to transform direction cosines of the
    #       particle back to the physical frame (from the transport frame)"

    # unlike original Mortran calculate phi directly.
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
