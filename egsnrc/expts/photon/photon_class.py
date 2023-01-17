# Photon class
"""Simple experiment to see if can pass class instance to GPU device functions"""

from dataclasses import dataclass
import numpy as np
import logging
import numba as nb

numba_logger = logging.getLogger('numba')
numba_logger.setLevel(logging.WARNING)

from numba import cuda


THREADS_PER_BLOCK = 10
BLOCKS = 2

@dataclass
class Particle:
    status: np.int32
    x: np.float32
    y: np.float32
    z: np.float32

@cuda.jit(device=True)
def update_particle(p, i):
    p.status = i


# Kernel
@cuda.jit
def test_kernel(particles, out):
    """Simple test to see what can be passed to numba cuda device function"""
    i = cuda.grid(1)
    if i < len(particles):
        p = Particle(*particles[i])
        update_particle(p, i)
        out[i] = p.status, p.x, p.y, p.z


if __name__ == "__main__":
    import sys
    NUM_PHOTONS = BLOCKS * THREADS_PER_BLOCK
    major, minor = sys.version_info.major, sys.version_info.minor
    print(f"Starting run with Numba {nb.__version__}, Python {major}.{minor}")
    print(f"Running {NUM_PHOTONS} photons")
    print(f"{BLOCKS=}   {THREADS_PER_BLOCK=}")
    # energies = [random_f32(rng_states, ) for i in range(100_000)]
    # energies = cuda.device_array(np.ones(1000000, dtype=np.float32))
    # energies_np = np.ones(1_000_000, dtype=np.float32)

    particles = np.zeros((NUM_PHOTONS, 4))
    out = np.empty_like(particles)
    particles[:,3] = 1.0  # Set z values to 1 just for kicks

    cuda.synchronize()
    test_kernel[BLOCKS, THREADS_PER_BLOCK](particles, out)
    cuda.synchronize()
    print(out)
