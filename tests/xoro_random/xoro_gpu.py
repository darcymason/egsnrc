# xoro_gpu.py
"""Test random number generation
"""
import numpy as np
import sys
from egsnrc import config
if int(sys.argv[-1]) == 64:
    config.KFLOAT = np.float64
    config.KINT = np.int64
else:
    config.KFLOAT = np.float32
    config.KINT = np.int32


from time import perf_counter
import numba as nb
from numba import cuda
import numpy as np
import sys
from egsnrc.config import KFLOAT, KINT
from egsnrc import egsrandom
from egsnrc.util import CUDATimer, cuda_details

import logging

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)


@cuda.jit
def random_test(rng_states, num_particles, out):
    """Fill in the output array with random [0, 1) energies"""
    # Get unique grid index
    gid = cuda.grid(1)

    # Don't need this check here because implicit in grid_stride loop
    # if gid >= num_particles:
    #     return

    threads_per_grid = cuda.blockDim.x * cuda.gridDim.x
    # if gid==0:
    #     print("threads_per_grid=", threads_per_grid)

    # CUDA Grid-Stride loop to reuse this thread
    for i_arr in range(gid, num_particles, threads_per_grid):
        egsrandom.random_kfloat(rng_states, gid)

    # Use last random so loop is not optimized away ... (?)
    if cuda.threadIdx.x == 0:
        out[cuda.blockIdx.x] = egsrandom.random_kfloat(rng_states, gid)


seed = 1
def run(num_particles, blocks_per_grid, threads_per_block):
    if not cuda.is_available():
        print("***** This script requires CUDA gpu  ****")
        sys.exit(-1)

    egsrandom.set_array_library("cuda")
    rng_states = egsrandom.initialize(seed, num_blocks * threads_per_block)

    cuda.synchronize()
    out = cuda.device_array(num_blocks, dtype=KFLOAT)
    # out = cuda.device_array(10, dtype=KFLOAT)  # XXX just a dummy array for now

    # Do a short run to jit the function so not included in the timing
    print("\n=============================")
    start = perf_counter()
    random_test[1, 10](
        rng_states, 10, out
    )
    cuda.synchronize()
    end = perf_counter()

    print(f"Initial short run to jit the function...{(end - start):>8.5} seconds")
    bits = 32 if KFLOAT is np.float32 else 64

    print("\n-------------------------------------------------------------")
    print(cuda_details())
    print(f"{bits}-bits run with {num_particles=:,} {blocks_per_grid=}  {threads_per_block=}")
    start2 = perf_counter()
    with CUDATimer() as cudatimer:  # CUDATimes(stream)
        random_test[num_blocks, threads_per_block](
            rng_states, num_particles, out
        )
    # cuda.synchronize()
    cuda_elapsed = cudatimer.elapsed
    end2 = perf_counter()
    elapsed_time = end2 - start2
    print(
        f"Elapsed time in ms: CUDA events: {cuda_elapsed:.2f};   "
        f"Python perf_counter {elapsed_time * 1000:.1f}"
    )
    print("\n-------------------------------------------------------------")
    print("\n Sample of random numbers out")
    out = out.copy_to_host()
    print(out[:30])
    print(out.dtype)
    print("First number, to show the precision:", out[0])




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