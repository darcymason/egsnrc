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
from egsnrc.util import CUDATimer

import logging

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)


@cuda.jit
def random_test(rng_states, num_particles, out):
    """Fill in the output array with random [0, 1) energies"""
    # Get unique grid index
    gid = cuda.grid(1)

    if gid >= num_particles:
        return

    # CUDA Grid-Stride loop to reuse this thread
    threads_per_grid = cuda.blockDim.x * cuda.gridDim.x
    for i_arr in range(gid, num_particles, threads_per_grid):
        egsrandom.random_kfloat(rng_states, gid)

    # Use last random so loop is not optimized away ... (?)
    out[gid] = egsrandom.random_kfloat(rng_states, gid)


seed = 1
def run(num_particles, blocks_per_grid, threads_per_block):
    print("\n-------------------------------")
    print(f"Run with {num_particles=:,} {blocks_per_grid=}  {threads_per_block=}\n")
    if not cuda.is_available():
        print("***** This script requires CUDA gpu  ****")
        sys.exit(-1)

    egsrandom.set_array_library("cuda")
    rng_states = egsrandom.initialize(seed, num_blocks * threads_per_block)

    print("Without storing the random numbers...")
    cuda.synchronize()
    out = cuda.device_array(num_blocks * threads_per_block, dtype=KFLOAT)
    # out = cuda.device_array(10, dtype=KFLOAT)  # XXX just a dummy array for now

    # Do a short run to jit the function so not included in the timing
    print("Initial short run to jit the function...")
    start = perf_counter()
    random_test[1, threads_per_block](
        rng_states, num_particles, out
    )
    end = perf_counter()
    print(f"Jit time (Python perf_counter): {(end - start):>8.5} seconds")
    print("Done jit.\n")
    start = perf_counter()
    with CUDATimer() as cudatimer:  # CUDATimes(stream)
        random_test[num_blocks, threads_per_block](
            rng_states, num_particles, out
        )
    cuda.synchronize()
    print(f"Elapsed time by events {cudatimer.elapsed:.2f} ms")
    elapsed_time = perf_counter() - start
    print(f"Elapsed time (Python perf_counter): {elapsed_time:>8.5} seconds")


if __name__ == "__main__":
    import sys
    from math import ceil

    num_particles = 2**20
    num_blocks = 512
    threads_per_block = 512

    # Int args could be num_particles, num_blocks, threads_per_block

    int_args = [int(x) for x in sys.argv[1:]] if len(sys.argv) > 1 else None

    if int_args:
        num_particles = int_args[0]
    try:
        num_blocks = int_args[1]
        threads_per_block = int_args[2]
    except:
        pass

    run(num_particles, num_blocks, threads_per_block)