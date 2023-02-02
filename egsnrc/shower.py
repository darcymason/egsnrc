
from egsnrc.photon import photon_kernel
from egsnrc import photon  # To set index in CPU mode
from egsnrc.config import on_gpu
from egsnrc import random
from egsnrc.media import Vacuum


def cuda_details():
    try:
        from cuda.cuda import (
            CUdevice_attribute, cuDeviceGetAttribute, cuDeviceGetName, cuInit
        )
    except ImportError:
        return (
            "** GPU details not available\n"
            "In Colab, use `!pip install cuda-python` to see GPU specs\n"
        )

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

def shower(
    seed, num_particles, source, regions, media, howfar, ausgab, iscore, fscore
):
    # Initialize random numbers
    if on_gpu:
        random.set_array_library("cuda")
    else:
        random.set_array_library("numpy")
    rng_states = random.initialize(seed)

    # Convert regions from Python classes to tuples needed to pass to kernel

    iparticles, fparticles = source.generate(num_particles)

    # Add vaccuum as medium 0 if not already there
    # XXX currently assumes media and regions are in medium/region # consecutive order
    # insert vaccum as medium 0 if not already there
    if media[0].number != 0:
        media = (Vacuum, *media)

    if not isinstance(media, tuple):  # necessary for current Numba kernels
        media = tuple(media)

    if regions[0].number != 0:
        regions = (tuple(), *regions)

    if not isinstance(regions, tuple):
        regions = tuple(regions)

    if on_gpu:
        photon_kernel.forall(num_particles)(
            rng_states,
            iparticles, fparticles,
            regions, media, howfar, ausgab, iscore, fscore
        )
    else:
        for i in range(num_particles):
            photon.non_gpu_index = i
            photon_kernel.py_func(
                rng_states,
                iparticles, fparticles,
                regions, media, howfar, ausgab, iscore, fscore
            )