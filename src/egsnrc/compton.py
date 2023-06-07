from math import log, exp, sqrt, sin, cos
from egsnrc.egsrandom import random_kfloat
from egsnrc.angles import uphi
from egsnrc.constants import REST_MASS, TWO_PI
from egsnrc.particles import replace_e_uvw
from egsnrc import config
from numba import cuda


@config.device_jit
def compton(rng_states, gid, p):
    """Return a modified particle from Compton interaction

    Returns
    -------
    mod_p   Particle
        The modified particle (new energy, direction cosines)
    """
    # TODO: do not use bound, always K-N
    # TODO: no e- created yet

    ko = p.energy / REST_MASS  # Gamma energy in units of electron rest energy
    broi = 1 + 2 * ko  # Needed for scattering angle sampling

    if ko > 2:  # At high energies the original EGS4 method is most efficient
        broi2 = broi * broi
        alph1 = log(broi)
        bro   = 1 / broi
        # alph2 = ko * (broi+1) * bro * bro # not used again, just inline below
        alpha = alph1 + ko * (broi+1) * bro * bro

        # Set up fake True for first pass through loop
        rnno19 = aux = br = 1
        rejf3 = 0.0
        while rnno19 * aux > rejf3 or not (bro < br < 1):  # rejection sampling loop
            rnno15 = random_kfloat(rng_states, gid)
            rnno16 = random_kfloat(rng_states, gid)
            rnno19 = random_kfloat(rng_states, gid)
            if rnno15 * alpha < alph1:  # Use 1/br part
                br = exp(alph1 * rnno16) * bro
            else:  # Use the br part
                br = sqrt(rnno16 * broi2 + (1 - rnno16))  * bro

            temp = (1 - br) / (ko * br)
            sinthe = max(0., temp*(2-temp))
            aux = 1 + br * br
            rejf3 = aux - br * sinthe

            # IF( br < 0.99999/broi | br > 1.00001 )
            #     $egs_warning(*,' sampled br outside of allowed range! ',ko,1./broi,br)

    else:  # At low energies it is faster to sample br uniformly
        bro = 1. / broi
        bro1 = 1 - bro
        rejmax = broi + bro

        # Set up fake True for first pass through loop
        rnno16 = br = 1.0
        rejf3 = 0.0
        while rnno16 * br * rejmax > rejf3 or not (bro < br < 1):
            rnno15 = random_kfloat(rng_states, gid)
            rnno16 = random_kfloat(rng_states, gid)
            br = bro + bro1 * rnno15
            temp = (1 - br) / (ko * br)
            sinthe = max(0., temp*(2-temp))
            rejf3 = 1 + br * br - br * sinthe
            # IF( br < 0.99999/broi | br > 1.00001 )
            #     $egs_warning(*,' sampled br outside of allowed range! ',ko,1./broi,br)

    # $RADC_REJECTION

    costhe = 1 - temp
    sinthe = sqrt(sinthe)
    energy = p.energy * br

    # Random sample the azimuth
    u, v, w = uphi(rng_states, gid, p, sinthe, costhe)

    # aux = 1 + br*br - 2*br*costhe
    # if aux > 1e-8:
    #     costhe = (1-br*costhe)/sqrt(aux)
    #     sinthe = (1-costhe)*(1+costhe)
    #     sinthe = -sqrt(sinthe) if sinthe > 0 else 0
    # else:
    #     costhe = 0
    #     sinthe = -1

    return replace_e_uvw(p, energy, u, v, w)
