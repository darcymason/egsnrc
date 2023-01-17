from math import log, exp, sqrt, sin, cos
import numpy as np
from numba import cuda
import numba as nb

from numba.cuda.random import create_xoroshiro128p_states
from numba.cuda.random import xoroshiro128p_uniform_float32 as random_f32
from collections import namedtuple

import logging


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


# Try to bundle attributes for easier use
# Note: in future should be able to combine all in one numba jitclass
#      (currently available only on CPU)

from numba import int32
# from typing import NamedTuple

# class iParticle(NamedTuple):
#     status: int32
#     region: int32


STATUS, REGION = range(2)
ENERGY, Z = range(2)

# Current particle Interaction
STEPPING, COMPTON, PHOTO = range(3)
SCORE_nCOMPTON, SCORE_nPHOTO, SCORE_eCOMPTON, SCORE_ePHOTO, SCORE_nLost, SCORE_eLost = range(6)

particle_iattrs = "status region".split()
particle_fattrs = "energy z".split()
iParticle = namedtuple("iParticle", particle_iattrs)
# iParticleType = NamedTuple((int32, int32), iParticle)

fParticle = namedtuple("fParticle", particle_fattrs)

# XXX dummy interaction probabilities
# compton_prob = cuda.to_device(np.array([0.8, 0.5], dtype=np.float32))

@cuda.jit(device=True)  # "Device function"
def scoring(gid, status1, region1, energy1, status2, region2, energy2, out):
    # Note, if don't track by gid, then need "atomic" operations to handle multi-thread
    if status1 == COMPTON:
        out[gid, region2, SCORE_nCOMPTON] += 1.0
        out[gid, region2, SCORE_eCOMPTON] += energy1 - energy2
    elif status1 == PHOTO:
        out[gid, region2, SCORE_nPHOTO] += 1.0
        out[gid, region2, SCORE_ePHOTO] += energy1 - energy2
    elif region2 == NUM_REGIONS:  # lost to geometry (0-based region index)
        out[gid, region2, SCORE_nLost] += 1.0
        out[gid, region2, SCORE_eLost] += energy1


@cuda.jit  # Kernel
def particle_kernel(rng_states, iparticles, fparticles, out):
    """Main particle simulation loop"""

    if not DEBUGGING_ON_CPU:
        gid = cuda.grid(1)
    if gid > len(fparticles):   # needed when have more GPU threads than particles
        return

    # oi = iparticles[gid]
    # of = fparticles[gid]
    # pi = iParticle(int32(oi[0]), int32(oi[1]))
    # pf = fParticle(of[0], of[1])

    status = iparticles[gid, STATUS]
    region = iparticles[gid, REGION]
    energy = fparticles[gid, ENERGY]
    z = fparticles[gid, Z]

    while status >= 0:  # Negative numbers for particle no longer tracked
        # Take a step
        status2 = STEPPING
        distance = random_f32(rng_states, gid)

        # Call "howfar"
        z2 = z + distance
        region2 = 1 if z2 >= 2.0 else 0
        if z2 >= 4:
            region2 = NUM_REGIONS # outside geometry (0-based counting)
            status2 = -1  # particle has left the geometry

        if status2 >= 0:
            rand = random_f32(rng_states, gid)
            if (region2 == 0 and rand < 0.8) or (region2 == 1 and rand < 0.5):
                # "Compton"
                status = COMPTON
                energy2 = energy * 0.8
                if energy2 < 0.001:  # PCUT
                    energy2 = 0
                    status2 = -2
            else:  # "photoelectric"
                status = PHOTO
                energy2 = 0
                status2 = -3

        # pi2 = iParticle(status, region)
        # pf2 = fParticle(energy, z)
        scoring(gid, status, region, energy, status2, region2, energy2, out)
        status, region = status2, region2
        energy, z = energy2, z2


if __name__ == "__main__":

    # To Debug on CPU:
    #  - set the flag True below
    # Add `gid` parameter to front of kernel call
    # Comment out the jit decorator for scoring function
    DEBUGGING_ON_CPU = False  # True
    from time import perf_counter
    import sys

    # THREADS_PER_BLOCK = 5
    # BLOCKS = 2
    if len(sys.argv) > 1:
        NUM_PHOTONS = int(sys.argv[1])
    NUM_REGIONS = 2

    Py_major, Py_minor = sys.version_info.major, sys.version_info.minor
    print(f"Starting run with Numba {nb.__version__}, Python {Py_major}.{Py_minor}")
    print(f"Running {NUM_PHOTONS:,} particles")
    print(cuda_details())
    # print(f"{BLOCKS=}   {THREADS_PER_BLOCK=}")

    if not DEBUGGING_ON_CPU:
        rng_states = create_xoroshiro128p_states(NUM_PHOTONS, seed=1)
        dev_rng_states = cuda.to_device(rng_states)
    times = []
    for run in range(3):
        energies_np = np.random.random(NUM_PHOTONS).astype(np.float32) + 0.511 # catch both ko>2 and <2
        fparticles = np.zeros((NUM_PHOTONS, len(particle_fattrs)))
        iparticles = np.zeros((NUM_PHOTONS, len(particle_iattrs)), dtype=np.int32)
        fparticles[:, ENERGY] = energies_np
        out = np.zeros(
            # 3 regions: 0, 1, and  2=escaped the geometry (just use eCompt)
            (NUM_PHOTONS, NUM_REGIONS + 1,  6),  # 6 for (nCompt, nPhoto, eCompt, ePhoto, nLost, eLost)
            dtype=np.float32
        )
        if not DEBUGGING_ON_CPU:
            dev_fparticles = cuda.to_device(fparticles)
            dev_iparticles = cuda.to_device(iparticles)
            dev_out = cuda.to_device(out)

        # Run
        start = perf_counter()
        if not DEBUGGING_ON_CPU:
            particle_kernel.forall(len(fparticles))(
                dev_rng_states, dev_iparticles, dev_fparticles, dev_out
            )
        else:
            random_f32 = lambda r,i: np.random.random(1)
            for i in range(NUM_PHOTONS):
                particle_kernel.py_func(i, None, iparticles, fparticles, out)
        end = perf_counter()

        times.append(end - start)

    print("----------------------")
    print("Times:", ', '.join(f"{time_:>8.5} " for time_ in times), "seconds")
    if not DEBUGGING_ON_CPU:
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
