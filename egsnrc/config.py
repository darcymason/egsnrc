from numba import cuda
from logging import getLogger

logger = getLogger("egsnrc")


def cuda_device_jit(func):
    return cuda.jit(func, device=True)


def no_jit(func):
    return func


def use_gpu(want_gpu=True):
    global device_jit
    global on_gpu

    device_jit = no_jit
    on_gpu = False
    if not want_gpu:
        return

    if cuda.is_available():
        device_jit = cuda_device_jit
        on_gpu = True
    else:
        logger.info("GPU not available, continuing with CPU only")


# Default is to use GPU if available, user can call use_gpu(False) if needed
use_gpu()
