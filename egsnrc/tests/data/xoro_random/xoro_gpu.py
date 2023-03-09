# example1.py
"""Define a slab geometry with photon tracking only
"""
# import os
# os.environ['NUMBA_ENABLE_CUDASIM'] = '1'  # XXX temp for testing GPU running on CPU

from time import perf_counter
import numba as nb
from numba import cuda
import numpy as np
import sys
from egsnrc.config import KFLOAT, KINT
from egsnrc import egsrandom
from egsnrc.particles import uniform_energies
from egsnrc.util import CUDATimer

import logging

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)


seed = 1
def run(num_particles, num_batches):
    if not cuda.is_available():
        print("***** This script requires CUDA gpu  ****")
        sys.exit(-1)

    egsrandom.set_array_library("cuda")
    rng_states = egsrandom.initialize(seed, num_particles)

    print("Without storing the random numbers...")
    cuda.synchronize()
    # XXX out = cuda.device_array(num_particles, dtype=KFLOAT)
    # out = cuda.device_array(10, dtype=KFLOAT)  # XXX just a dummy array for now

    # Do a short run to jit the function so not included in the timing
    print("Initial short run to jit the function...")
    start = perf_counter()
    uniform_energies.forall(KINT(100))(
        rng_states, num_particles, # out
    )
    end = perf_counter()
    print(f"Jit time (Python perf_counter): {(end - start):>8.5} seconds")
    print("Done jit.")
    start = perf_counter()
    with CUDATimer() as cudatimer:  # CUDATimes(stream)
        uniform_energies.forall(KINT(num_particles))(
            rng_states, num_particles #, out
        )
    cuda.synchronize()
    print(f"Elapsed time by events {cudatimer.elapsed:.2f} ms")
    elapsed_time = perf_counter() - start
    print(f"Elapsed time (Python perf_counter): {elapsed_time:>8.5} seconds")


if __name__ == "__main__":
    import sys
    num_particles = 200
    num_batches = 8
    if len(sys.argv) > 1:
        try:
            num_particles = int(sys.argv[1])
        except ValueError:
            print("Optional num_particles command-line argument must be an integer")
            print(usage)
            sys.exit(-1)
    if len(sys.argv) > 2:
        try:
            num_batches = int(sys.argv[2])
        except ValueError:
            print("Optional num_batches command-line argument must be an integer")
            print(usage)

    run(num_particles, num_batches)