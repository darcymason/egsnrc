from math import log, exp, sqrt, sin, cos
import numpy as np
from numba import cuda
import numba as nb

from numba.cuda.random import create_xoroshiro128p_states
from numba.cuda.random import xoroshiro128p_uniform_float32 as random_f32
from numba import int32, float32
from numba.core.types import NamedTuple

from collections import namedtuple


import logging


try:
    from cuda.cuda import CUdevice_attribute, cuDeviceGetAttribute, cuDeviceGetName, cuInit
    have_cuda_python = True
except ImportError:
    have_cuda_python = False

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)


particle_iattrs = "status region".split()
particle_fattrs = "energy z".split()


STATUS, REGION = np.arange(2, dtype=np.int32)
ENERGY, Z = np.arange(2, dtype=np.int32)

# Current particle Interaction
STEPPING, COMPTON, PHOTO = np.arange(3, dtype=np.int32)
SCORE_nCOMPTON, SCORE_nPHOTO, SCORE_eCOMPTON, SCORE_ePHOTO, \
    SCORE_nLost, SCORE_eLost = np.arange(6, dtype=np.int32)

Particle = namedtuple("Particle", ("status", "region", "energy", "z"))
ParticleType = NamedTuple((int32, int32, float32, float32), Particle)


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



# XXX dummy interaction probabilities
# compton_prob = cuda.to_device(np.array([0.8, 0.5], dtype=np.float32))

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
def particle_kernel(rng_states, iparticles, fparticles, p, p2, out):
    """Main particle simulation loop"""

    if not DEBUGGING_ON_CPU:
        gid = cuda.grid(1)  # grid index - unique in entire grid
    if gid > len(fparticles):   # needed when have more GPU threads than particles
        return

    # Pack array info into a Particle namedtuple, for convenience
    # Note tuples are not mutable, so e.g. `status`` updates must stand outside
    oi = iparticles[gid]
    of = fparticles[gid]
    p = Particle(oi[STATUS], oi[REGION], of[ENERGY], of[Z])

    status = p.status  # see note above

    while status >= 0:  # Negative numbers for particle no longer tracked
        # Take a step
        status2 = nb.int32(STEPPING)
        distance = random_f32(rng_states, gid)

        # Call "howfar"
        z2 = float32(p.z + distance)
        region2 = int32(1) if z2 >= 2.0 else int32(0)
        if z2 >= 4.0:
            region2 = NUM_REGIONS # outside geometry (0-based counting)
            status2 = int32(-1)  # particle has left the geometry

        if status2 >= 0:
            rand = random_f32(rng_states, gid)
            if (region2 == 0 and rand < 0.8) or (region2 == 1 and rand < 0.5):
                # "Compton"
                status = COMPTON
                energy2 = float32(p.energy * 0.8)
                if energy2 < 0.001:  # PCUT
                    energy2 = float32(0)
                    status2 = int32(-2)
            else:  # "photoelectric"
                status = PHOTO
                energy2 = float32(0)
                status2 = int32(-3)

        p2 = Particle(status2, region2, energy2, z2)
        # p = p._replace(status=status)  # replace doesn't work in Cuda
        p = Particle(status, p.region, p.energy, p.z)
        scoring(gid, p, p2, out)  # probs don't need gid, can get with cuda.grid(1)
        # New particle info becomes current for next loop
        p = p2
        status = p.status  # need mutable status


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
        num_photons = int(sys.argv[1])
    else:
        num_photons = 20  # for testing
    NUM_REGIONS = int32(2)

    Py_major, Py_minor = sys.version_info.major, sys.version_info.minor
    print(f"Starting run with Numba {nb.__version__}, Python {Py_major}.{Py_minor}")
    print(f"Running {num_photons:,} particles")
    print(cuda_details())
    if DEBUGGING_ON_CPU:
        print("**** NOTE: DEBUGGING_ON_CPU is on  ****")

    # print(f"{BLOCKS=}   {THREADS_PER_BLOCK=}")

    if not DEBUGGING_ON_CPU:
        rng_states = create_xoroshiro128p_states(num_photons, seed=1)
        dev_rng_states = cuda.to_device(rng_states)
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
        if not DEBUGGING_ON_CPU:
            dev_fparticles = cuda.to_device(fparticles)
            dev_iparticles = cuda.to_device(iparticles)
            dev_out = cuda.to_device(out)

        # Run
        if not DEBUGGING_ON_CPU:
            cuda.synchronize()
        start = perf_counter()
        if not DEBUGGING_ON_CPU:
            # Try to force typing
            p = Particle(int32(0), int32(0), float32(1), float32(0))
            p2 = p
            print(f"{type(iparticles[0, 0])=}")
            print(f"{type(fparticles[0, 0])=}")
            particle_kernel.forall(len(fparticles))(
                dev_rng_states, dev_iparticles, dev_fparticles, p, p2, dev_out
            )
        else:
            random_f32 = lambda r,i: np.random.random(1)
            for i in range(num_photons):
                particle_kernel.py_func(i, None, iparticles, fparticles, out)
        if not DEBUGGING_ON_CPU:
            cuda.synchronize()
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
