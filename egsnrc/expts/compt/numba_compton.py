from math import log, exp, sqrt, sin, cos
# from egsnrc import random
from egsnrc.constants import RM, TWO_PI
import numpy as np
from numba import cuda
import numba as nb
import random as pyrandom

import logging

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

@nb.jit(debug=False, nopython=True)
def numba_compton_jit(energy):
    """Subroutine for sampling incoherent (Compton) scattering"""
    # calc_azimith False used in test_compton only
    # TODO: do not use bound, always K-N
    # TODO: no e- created yet

    ko = energy / RM  # Gamma energy in units of electron rest energy
    broi = 1 + 2 * ko  # Needed for scattering angle sampling

    if ko > 2:  # At high energies the original EGS4 method is most efficient
        broi2 = broi * broi
        alph1 = log(broi)
        bro   = 1 / broi
        # alph2 = ko * (broi+1) * bro * bro # not used again, just inline below
        alpha = alph1 + ko * (broi+1) * bro * bro

        # Set up fake True for first pass through loop
        rnno19 = aux = br = 1; rejf3 = 0.0
        while rnno19 * aux > rejf3 or not (bro < br < 1):  # rejection sampling loop
            rnno15 = pyrandom.random()
            rnno16 = pyrandom.random()
            rnno19 = pyrandom.random()
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
        rnno16 = br = 1.0; rejf3 = 0.0
        while rnno16 * br * rejmax > rejf3 or not (bro < br < 1):
            rnno15 = pyrandom.random()
            rnno16 = pyrandom.random()
            br = bro + bro1 * rnno15
            temp = (1 - br) / (ko * br)
            sinthe = max(0., temp*(2-temp))
            rejf3 = 1 + br * br - br * sinthe
            # IF( br < 0.99999/broi | br > 1.00001 )
            #     $egs_warning(*,' sampled br outside of allowed range! ',ko,1./broi,br)

    # $RADC_REJECTION

    costhe = 1 - temp
    sinthe = sqrt(sinthe)
    energy *= br
    # Random sample the azimuth
    # if calc_azimuth:
    #     _, (azimuth_ran,) = random.floats_0_1(rng, 1)
    #     phi = TWO_PI * azimuth_ran
    #     sinphi = sin(phi)
    #     cosphi = cos(phi)
    # else:
    sinphi = cosphi = 99  # XXX temp
    # aux = 1 + br*br - 2*br*costhe
    # if aux > 1e-8:
    #     costhe = (1-br*costhe)/sqrt(aux)
    #     sinthe = (1-costhe)*(1+costhe)
    #     sinthe = -sqrt(sinthe) if sinthe > 0 else 0
    # else:
    #     costhe = 0
    #     sinthe = -1

    return energy, sinthe, costhe, sinphi, cosphi

i32 = nb.int32
# nb.jit("i32(i32)->i32", nopython=True)
def xxx_numba(n):
    return sum((i + i**2 + sqrt(i**2) for i in range(n)))


if __name__ == "__main__":
    from time import perf_counter
    # random.set_array_library("numpy")
    pyrandom.seed(42)
    energies = [pyrandom.random() for i in range(100_000)]
    energies = np.array(energies) + 0.001
    times = []
    for run in range(3):
        from time import perf_counter
        start = perf_counter()
        for energy in energies:
            numba_compton_jit(energy)
        end = perf_counter()
        times.append(end - start)
    print("----------------------")
    print("Times:", ', '.join(f"{time_:>8.5} " for time_ in times), "seconds")
