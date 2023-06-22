import numpy as np
import torch

from egsnrc.particles import Particles
from egsnrc import egsrandom


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

    egsrandom.set_array_library(array_library)
    key = egsrandom.initialize(42)

    particles = Particles(f_arr, i_arr, array_library=array_library)  # device?

    return particles, key, device

# ======================================
# PHOTON "transport"
def run(particles, key, device):
    # Start in region 0, medium 0, status 0 so leave those alone
    particles.f_arr[Particles.ENERGY] = 100.0  # particles.energy = 100.0
    particles.f_arr[Particles.W] = 1.0

    while particles.any_alive():
        particles.log_energy()
        # Sample number of mfp to transport before interacting
        rnno35 = egsrandom.float_0_1(len(particles), device)
        rnno35[rnno35==0] += 1.0e-30  # could use `numpy.clip` but not quite same as old code
        torch.log(rnno35)
        key, ran_floats = egsrandom.floats_0_1(key, len(particles.energy), device)
        particles.f_arr[Particles.ENERGY] -= 20.0 * ran_floats


def time_run(num_particles, array_lib, want_gpu, repeats=10):
    # Simple test
    import timeit
    print("Starting time_run...")
    particles, key, device = setup(num_particles, array_lib, want_gpu)
    print("Particles created")
    time_cuda = None
    if device.startswith("cuda"):
        print("Have cuda. Starting timer")
        torch.cuda.synchronize()  # wait for any operations to complete
        start = torch.cuda.Event(enable_timing=True)
        end = torch.cuda.Event(enable_timing=True)
        start.record()
    print(f"Starting runs with {array_lib=}  {device=}...")
    run_func = run
    seconds = timeit.timeit(
        "run_func(particles, key, device=device)", globals=locals(), number=repeats
    )
    if device.startswith("cuda"):
        print("Stopping cuda timer")
        torch.cuda.synchronize()  # wait for all_reduce to complete
        end.record()
        torch.cuda.synchronize()  # need to wait once more for op to finish
        time_cuda = start.elapsed_time(end) / 1000

    print(f"Timeit (repeats={repeats}): {seconds} seconds")
    print(f"{time_cuda=}")

    # run(particles, key, device=device)

if __name__ == "__main__":
    time_run(10_000, "pytorch", False)

