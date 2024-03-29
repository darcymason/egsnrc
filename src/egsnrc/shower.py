import sys
from time import perf_counter
from numba import cuda
import numba as nb
import logging

try:
    from cuda.cuda import (
        CUdevice_attribute,
        cuDeviceGetAttribute,
        cuDeviceGetName,
        cuInit,
    )

    have_cuda_python = True
except ImportError:
    have_cuda_python = False


from egsnrc.photon import photon_kernel
from egsnrc import photon  # To set index in CPU mode and howfar, ausgab
from egsnrc.config import on_gpu
from egsnrc import egsrandom
from egsnrc.media import Vacuum
from egsnrc.util import CUDATimer

logger = logging.getLogger("egsnrc")


def cuda_details():
    try:
        from cuda.cuda import (
            CUdevice_attribute,
            cuDeviceGetAttribute,
            cuDeviceGetName,
            cuInit,
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

    attrs = {"DEVICE_NAME": DEVICE_NAME.strip()}
    attr_names = [
        "MAX_THREADS_PER_BLOCK",
        "MAX_BLOCK_DIM_X",
        "MAX_GRID_DIM_X",
        "MULTIPROCESSOR_COUNT"
    ]
    for attr in attr_names:
        err, attrs[attr] = cuDeviceGetAttribute(
            getattr(CUdevice_attribute, f"CU_DEVICE_ATTRIBUTE_{attr}"), 0
        )

    return attrs


threads_per_block = 512  # 1024
# blocks_per_grid = 16*40


def shower(
    seed,
    num_particles,
    get_source_particle,
    regions,
    media,
    howfar,
    ausgab,
    iscore,
    fscore,
):
    """Start the Monte Carlo simulation

    Parameters
    ----------
    seed : int
        Random number generator seed
    num_particles : int
        Number of photons to launch
    source : Source subclass
        A radition Source class with `generate` method
    regions : tuple
        Tuple of _Region namedtuples
    media : tuple
        Tuple of _Medium namedtuples
    iscore : array of num_particles int32's
        Output array for ausgab to track integer values, e.g. counting interactions
    fscore : array of num_particles float32's
        Output array for ausgab to track float values, e.g. energy deposited
    """
    # Initialize random numbers
    if on_gpu:
        egsrandom.set_array_library("cuda")
    else:
        egsrandom.set_array_library("numpy")
    rng_states = egsrandom.initialize(seed, num_particles)

    # iparticles, fparticles = source.generate(num_particles)
    device = "cuda" if on_gpu else "cpu"
    if device == "cuda":
        print("Pre-copying arrays to device - not included in timing")
        # iparticles = cuda.to_device(iparticles)
        # fparticles = cuda.to_device(fparticles)
        iscore = cuda.to_device(iscore)
        fscore = cuda.to_device(fscore)
        rng_states = cuda.to_device(rng_states)
        blocks_per_grid = (num_particles // threads_per_block) + 1
        print(f"{threads_per_block=}  {blocks_per_grid=}")
    # Convert regions from Python classes to tuples needed to pass to kernel
    # Add vaccuum as medium 0 if not already there
    # XXX currently assumes media and regions are in medium/region # consecutive order
    # XXX also duplicates memory via Media kernelize - in Media and in Region.medium
    # insert vacuum as medium 0 if not already there

    # Convert Media types to be kernel compatible
    media = [medium.kernelize() for medium in media]
    if media[0].number != 0:
        media = (Vacuum, *media)

    if not isinstance(media, tuple):  # necessary for current Numba kernels
        media = tuple(media)

    # Convert Region types to be kernel compatible
    regions = [region.kernelize() for region in regions]
    if regions[0].number != 0:
        regions = (tuple(), *regions)  # XXX need a proper Region here

    if not isinstance(regions, tuple):
        regions = tuple(regions)

    # Configure callbacks:
    photon.howfar = howfar
    photon.ausgab = ausgab
    photon.get_source_particle = get_source_particle
    Py_major, Py_minor = sys.version_info.major, sys.version_info.minor
    logger.info(
        f"Starting `shower` with Numba {nb.__version__}, Python {Py_major}.{Py_minor}"
    )
    logger.info(f"Running {num_particles:,} particles")
    logger.info(cuda_details())

    if on_gpu:
        cuda.synchronize()
        start = perf_counter()
        with CUDATimer() as cudatimer:  # CUDATimes(stream)
            photon_kernel[blocks_per_grid, threads_per_block](
                rng_states,
                num_particles,  # XXX temp
                # iparticles, fparticles,
                regions,
                media,
                iscore,
                fscore,
            )
        print(f"Elapsed time by events {cudatimer.elapsed:.2f} ms")
        cuda.synchronize()
        end = perf_counter()
    else:
        start = perf_counter()
        for i in range(num_particles):
            photon.non_gpu_index = i
            photon_kernel.py_func(
                rng_states,
                num_particles,
                # iparticles, fparticles,
                regions,
                ausgab,
                iscore,
                fscore,
            )
        end = perf_counter()

    shower_time = end - start

    logger.info(f"Elapsed time (Python perf_counter): {shower_time:>8.5} seconds")

    if device == "cuda":
        iscore = iscore.copy_to_host()
        fscore = fscore.copy_to_host()
    return iscore, fscore
