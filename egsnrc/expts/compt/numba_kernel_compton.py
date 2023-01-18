from math import log, exp, sqrt, sin, cos
# from egsnrc import random
from egsnrc.constants import RM, TWO_PI
import numpy as np
from numba import cuda
import numba as nb
import random as pyrandom

# from numba.cuda.random import create_xoroshiro128p_states
# from numba.cuda.random import xoroshiro128p_uniform_float32 as random_f32

import logging

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

THREADS_PER_BLOCK = 512
BLOCKS = 4096
# rng_states = create_xoroshiro128p_states(THREADS_PER_BLOCK * BLOCKS, seed=1)


def cuda_details():
    # Need to: pip install --upgrade cuda-python

    from cuda.cuda import CUdevice_attribute, cuDeviceGetAttribute, cuDeviceGetName, cuInit

    # Initialize CUDA Driver API
    (err,) = cuInit(0)

    # Get attributes
    err, DEVICE_NAME = cuDeviceGetName(128, 0)
    DEVICE_NAME = DEVICE_NAME.decode("ascii").replace("\x00", "")

    attrs = {'DEVICE_NAME': DEVICE_NAME.strip()}
    attr_names = "MAX_THREADS_PER_BLOCK MAX_BLOCK_DIM_X MAX_GRID_DIM_X MULTIPROCESSOR_COUNT".split()
    for attr in attr_names:
        err, attrs[attr] =  cuDeviceGetAttribute(
        getattr(CUdevice_attribute, f"CU_DEVICE_ATTRIBUTE_{attr}"), 0
    )

    return attrs


@cuda.jit
def numba_kernel_compton(rng_states, energy, out_energy, out_costhe):
    """Subroutine for sampling incoherent (Compton) scattering"""
    # calc_azimith False used in test_compton only
    # TODO: do not use bound, always K-N
    # TODO: no e- created yet
    i = cuda.grid(1)
    if i < len(energy):
        ko = energy[i] / RM  # Gamma energy in units of electron rest energy
        broi = 1 + 2 * ko  # Needed for scattering angle sampling

        if ko > 2:  # At high energies the original EGS4 method is most efficient
            broi2 = broi * broi
            alph1 = log(broi)
            bro   = 1 / broi
            # alph2 = ko * (broi+1) * bro * bro # not used again, just inline below
            alpha = alph1 + ko * (broi+1) * bro * bro

            # Set up fake True for first pass through loop
            rnno19 = aux = br = 1.0; rejf3 = 0.0
            iteration = 0
            while rnno19 * aux > rejf3 or not (bro < br < 1):  # rejection sampling loop
                # rnno15 = random_f32(rng_states, i+iteration)
                rnno15 = rng_states[i+iteration]
                rnno16 = rng_states[i+iteration+1]
                rnno19 = rng_states[i+iteration+2]
                iteration += 5  # Just something to move the random state
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
            iteration = 0
            while rnno16 * br * rejmax > rejf3 or not (bro < br < 1):
                # rnno15 = random_f32(rng_states, thread_id)
                rnno15 = rng_states[i+iteration]
                rnno16 = rng_states[i+iteration+1]
                iteration += 5
                br = bro + bro1 * rnno15
                temp = (1 - br) / (ko * br)
                sinthe = max(0., temp*(2-temp))
                rejf3 = 1 + br * br - br * sinthe
                # IF( br < 0.99999/broi | br > 1.00001 )
                #     $egs_warning(*,' sampled br outside of allowed range! ',ko,1./broi,br)

        # $RADC_REJECTION

        out_costhe[i] = 1 - temp
        sinthe = sqrt(sinthe)
        out_energy[i] = energy[i] * br
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

    #  return energy, sinthe, costhe, sinphi, cosphi # no return for kernel

i32 = nb.int32
# nb.jit("i32(i32)->i32", nopython=True)
def xxx_numba(n):
    return sum((i + i**2 + sqrt(i**2) for i in range(n)))


if __name__ == "__main__":
    from time import perf_counter
    import sys

    NUM_PHOTONS = BLOCKS * THREADS_PER_BLOCK
    major, minor = sys.version_info.major, sys.version_info.minor
    print(f"Starting run with Numba {nb.__version__}, Python {major}.{minor}")
    print(f"Running {NUM_PHOTONS} photons")
    print(cuda_details())
    print(f"{BLOCKS=}   {THREADS_PER_BLOCK=}")
    # energies = [random_f32(rng_states, ) for i in range(100_000)]
    # energies = cuda.device_array(np.ones(1000000, dtype=np.float32))
    # energies_np = np.ones(1_000_000, dtype=np.float32)

    times = []
    rng_states_np = np.random.random(NUM_PHOTONS + 500).astype(np.float32)
    dev_rng_states = cuda.to_device(rng_states_np)
    for run in range(3):
        energies_np = np.random.random(NUM_PHOTONS).astype(np.float32) + 0.511 # catch both ko>2 and <2
        dev_energies = cuda.to_device(energies_np)
        out_energy = cuda.device_array_like(dev_energies)
        out_costhe = cuda.device_array_like(dev_energies)
        cuda.synchronize()
        start = perf_counter()
        numba_kernel_compton[BLOCKS, THREADS_PER_BLOCK](dev_rng_states, dev_energies, out_energy, out_costhe)
        cuda.synchronize()
        end = perf_counter()
        times.append(end - start)
    print("----------------------")
    print("Times:", ', '.join(f"{time_:>8.5} " for time_ in times), "seconds")
    host_energies = out_energy.copy_to_host()
    host_costhe = out_costhe.copy_to_host()
    sample = 100
    print(f"First {sample} of various arrays:")
    print(f"  Input energies: {energies_np[:sample]}")
    print(f"  Input randoms : {rng_states_np[:sample]}")
    print(f"  Energy out    : {host_energies[:sample]}")
    print(f"  Costhe out    : {host_costhe[:sample]}")
