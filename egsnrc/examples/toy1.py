import numpy as np
import torch

from egsnrc.particles import Particles
from egsnrc import random


def setup(num_particles, array_library, want_gpu=False):
    if array_library == "numpy" and want_gpu:
        raise ValueError("No GPU available for numpy")

    if array_library == "pytorch":
        device = "cuda" if want_gpu else "cpu"
        if want_gpu and not torch.cuda.is_available():
            raise ValueError("GPU is not currently available")
    else:
        device = None

    P = Particles  # short-form for accessing indices
    # Set up toy slab geometry like tutor examples
    if array_library == "numpy":
        f_arr = np.zeros((P.NUM_FLOAT_PARAMS, num_particles) , dtype=np.float32)
        i_arr = np.zeros((P.NUM_INT_PARAMS, num_particles), dtype=np.int32)

    elif array_library == "pytorch":
        f_arr = torch.zeros(
            (P.NUM_FLOAT_PARAMS, num_particles) , dtype=torch.float32, device=device
        )
        i_arr = torch.zeros(
            (P.NUM_INT_PARAMS, num_particles), dtype=torch.int32, device=device
        )
    else:
        raise NotImplementedError("Unknown array library")

    random.set_array_library(array_library)
    key = random.initialize(42)

    particles = Particles(f_arr, i_arr, array_library=array_library)  # device?

    return particles, key, device

# ======================================
# PHOTON "transport"
def run(particles, key, device):
    # Start in region 0, medium 0, status 0 so leave those alone
    P = Particles

    particles.f_arr[P.ENERGY] = 100.0
    # f_arr[P.W, :] = 1.0

    while particles.any_alive():
        key, ran_floats = random.floats_0_1(key, len(particles.energy), device=device)
        particles.energy -= 20.0 * ran_floats


if __name__ == "__main__":
    # Simple test
    import timeit
    num_particles = 10_000
    # array_lib = "numpy"
    array_lib = "pytorch"
    want_gpu = True
    particles, key, device = setup(num_particles, array_lib, want_gpu)
    print(f"Starting runs with {array_lib=}  {device=}...")
    seconds = timeit.timeit(
        "run(particles, key, device=device)", globals=globals(), number=10
    )
    print(f"Completed run in {seconds} seconds")



    # run(particles, key, device=device)


