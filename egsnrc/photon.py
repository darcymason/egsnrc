from math import log, exp, sqrt, sin, cos
import sys
import numpy as np
from numba import cuda
import numba as nb

from numba.cuda.random import create_xoroshiro128p_states
from numba.cuda.random import xoroshiro128p_uniform_float32 as random_f32
from numba import int32, float32
from numba.core.types import NamedTuple

import logging

from egsnrc.particles import Particle, STATUS, REGION, ENERGY, Z
from egsnrc.particles import STEPPING, COMPTON, PHOTO
from egsnrc.particles import GEOMETRY_DISCARD, PCUT_DISCARD, PHOTO_DISCARD

try:
    from cuda.cuda import CUdevice_attribute, cuDeviceGetAttribute, cuDeviceGetName, cuInit
    have_cuda_python = True
except ImportError:
    have_cuda_python = False

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)



def cuda_details():
    if not have_cuda_python:
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



# Example of scoring function (ausgab) - replace with user one
@cuda.jit(device=True)  # "Device function"
def scoring(gid, p1, p2, out):
    # Note, if don't track by gid, then need "atomic" operations to handle multi-thread
    if p1.status == COMPTON:
        out[gid, p2.region, SCORE_nCOMPTON] += 1.0
        out[gid, p2.region, SCORE_eCOMPTON] += p1.energy - p2.energy
    elif p1.status == PHOTO:
        out[gid, p2.region, SCORE_nPHOTO] += 1.0
        out[gid, p2.region, SCORE_ePHOTO] += p1.energy - p2.energy
    elif p2.region == NUM_REGIONS:  # lost to geometry (0-based region index)
        out[gid, p2.region, SCORE_nLost] += 1.0
        out[gid, p2.region, SCORE_eLost] += p1.energy


# Kernel
@cuda.jit
def photon_kernel(rng_states, iparticles, fparticles, regions, howfar, ausgab, out):
    """Main photon particle simulation kernel - implicitly called on each gpu thread"""

    gid = cuda.grid(1)  # grid index - unique in entire grid
    if gid > len(fparticles):   # needed when have more GPU threads than particles
        return

    # Pack array info into a Particle namedtuple, for convenience
    # Note tuples are not mutable, so e.g. `status`` updates must stand outside
    oi = iparticles[gid]
    of = fparticles[gid]
    p = Particle(oi[STATUS], oi[REGION], of[ENERGY], of[Z])

    status = p.status  # see note above, need mutable
    while status >= 0:  # Negative numbers for particle no longer tracked
        # Take a step
        status2 = nb.int32(STEPPING)
        distance = random_f32(rng_states, gid)

        region, discard = howfar(p, regions)

        # Are interacting
        if status2 >= 0:
            rand = random_f32(rng_states, gid)
            if (region2 == 0 and rand < 0.8) or (region2 == 1 and rand < 0.5):
                # "Compton"
                status = COMPTON
                energy2 = float32(p.energy * 0.8)
                if energy2 < 0.001:  # PCUT
                    energy2 = float32(0)
                    status2 = PCUT_DISCARD
            else:  # "photoelectric"
                status = PHOTO_DISCARD
                energy2 = float32(0)
                status2 = int32(-3)

        p2 = Particle(status2, region2, energy2, z2)
        # p = p._replace(status=status)  # replace doesn't work in Cuda
        p = Particle(status, p.region, p.energy, p.z)
        ausgab(gid, p, p2, out)  # probs don't need gid, can get with cuda.grid(1)
        # New particle info becomes current for next loop
        p = p2
        status = p.status  # need mutable status

def init(random_seed, num_particles):
    rng_states = create_xoroshiro128p_states(
        num_particles,
        seed=random_seed
    )
    return cuda.to_device(rng_states)


def run(particles, scoring_out, on_gpu=True):
    Py_major, Py_minor = sys.version_info.major, sys.version_info.minor
    print(f"Starting run with Numba {nb.__version__}, Python {Py_major}.{Py_minor}")
    print(f"Running {len(particles):,} particles")

    if cuda.is_available():
        print(cuda_details())
    else:
        thread_index = GridIterator()
        print("**** Running on CPU  ****")

    # To Debug on CPU:
    #  - set the flag True below
    # Add `gid` parameter to front of kernel call
    # Comment out the jit decorator for scoring function
    debug_on_cpu = False  # True
    from time import perf_counter


    particle_kernel.forall(len(fparticles))(
        dev_rng_states, dev_iparticles, dev_fparticles, dev_out
    )

    times = []
    for run in range(3):
        energies_np = np.random.random(num_photons) + 0.511 # catch both ko>2 and <2
        energies_np = energies_np.astype(np.float32)
        fparticles = np.zeros((num_photons, len(particle_fattrs)), dtype=np.float32)
        iparticles = np.zeros((num_photons, len(particle_iattrs)), dtype=np.int32)
        fparticles[:, ENERGY] = energies_np
        out = np.zeros(
            # 3 regions: 0, 1, and  2=escaped the geometry (just use eCompt)
            (num_photons, NUM_REGIONS + 1,  6),  # 6 for (nCompt, nPhoto, eCompt, ePhoto, nLost, eLost)
            dtype=np.float32
        )
        if not debug_on_cpu:
            dev_fparticles = cuda.to_device(fparticles)
            dev_iparticles = cuda.to_device(iparticles)
            dev_out = cuda.to_device(out)

        # Run
        if not debug_on_cpu:
            cuda.synchronize()
        start = perf_counter()
        if not debug_on_cpu:
            pass
            # Try to force typing
            # p = Particle(int32(0), int32(0), float32(1), float32(0))
            # p2 = p
            # print(f"{type(iparticles[0, 0])=}")
            # print(f"{type(fparticles[0, 0])=}")

        else:
            random_f32 = lambda r,i: np.random.random(1)
            for i in range(num_photons):
                particle_kernel.py_func(i, None, iparticles, fparticles, out)
        if not debug_on_cpu:
            cuda.synchronize()
        end = perf_counter()

        times.append(end - start)

    print("----------------------")
    print("Times:", ', '.join(f"{time_:>8.5} " for time_ in times), "seconds")
    if not debug_on_cpu:
        out = dev_out.copy_to_host()
    if len(out) < 50:
        print("Scoring - Particles/Region")
        print("nCompt     nPhoto     eCompt     ePhoto     nLost      eLost")
        print(out)
    for r in range(2):
        print("Energy Total by Region")
        print(f"Region {r}----")
        print("Compton: ", sum(out[:,r,SCORE_eCOMPTON]))
        print("Photo  : ", sum(out[:,r,SCORE_ePHOTO]))
    sum_compt = np.sum(out[:, :, SCORE_eCOMPTON])
    sum_photo = np.sum(out[:, :, SCORE_ePHOTO])
    sum_lost = np.sum(out[:, :, SCORE_eLost])
    print("Sum Compt, Photo, Lost", sum_compt, sum_photo, sum_lost)
    print("Energy in :", sum(fparticles[:, ENERGY]))
    print("Energy out:", sum((sum_compt, sum_photo, sum_lost)))

    # print("In particles: ---------------")
    # print(particles)
    # print("Out particles: ---------------")
    # print(out_particles)
    # sample = 100
    # print(f"First {sample} of various arrays:")
    # print(f"  Input energies: {energies_np[:sample]}")
    # print(f"  Input randoms : {rng_states_np[:sample]}")
    # print(f"  Energy out    : {host_energies[:sample]}")
    # print(f"  Costhe out    : {host_costhe[:sample]}")
