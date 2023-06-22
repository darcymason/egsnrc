import os
from numba import cuda
import numpy as np
from logging import getLogger

logger = getLogger("egsnrc")

try:
    cudasim_env = os.environ['NUMBA_ENABLE_CUDASIM']
except KeyError:
    cudasim_env = 0
on_cuda_sim = bool(int(cudasim_env))


# Define bit-width for all int's and floats in the GPU kernels
# NOTE: need to set egsrandom routines to same bits
KINT = np.int32
KFLOAT = np.float32


def cuda_device_jit(func):
    return cuda.jit(func, device=True)


def no_jit(func):
    return func


def use_gpu(want_gpu=True):
    global device_jit
    global on_gpu

    device_jit = no_jit  # nb.njit   # no_jit
    on_gpu = False
    if not want_gpu:
        return

    if cuda.is_available() or on_cuda_sim:
        device_jit = cuda_device_jit
        on_gpu = True
        if on_cuda_sim:
            logger.info("** Numba CUDA Simulator is in use **")

    else:
        logger.info("GPU not available, continuing with CPU only")


# Default is to use GPU if available, user can call use_gpu(False) if needed
use_gpu()
