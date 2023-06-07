
import itertools
import numpy as np
import numba

# NOTE: if change from 32-bit floats, need to change xoroshiro here
from numba.cuda.random import create_xoroshiro128p_states
from numba.cuda.random import xoroshiro128p_uniform_float32


def _np_initialize(seed, num_particles=None):
    return np.random.default_rng(seed)

def _np_float32(rng, _):
    return rng.random()

# Set up "random" from known sequence -------------
#    (e.g. to compare with mortran EGSnrc)
def _seq_initialize(known_list):
    return SeqGen(known_list)

def _seq_float32(rng, num=1, device=None):
    return rng.random()

def _cuda_initialize(seed, num_particles):
    return numba.cuda.to_device(create_xoroshiro128p_states(num_particles, seed))

class SeqGen:
    def __init__(self, known_list):
        self.known_list = known_list
        self.pos = 0

    def random(self):
        rand = self.known_list[self.pos]
        self.pos += 1
        return rand


# Configuration to select between libraries
def set_array_library(lib: str):
    global random_kfloat
    global initialize

    if lib == "numpy":
        random_kfloat = _np_float32  # Later switch depending on using gpu arrays
        initialize = _np_initialize
    elif lib == "sequence":
        random_kfloat = _seq_float32
        initialize = _seq_initialize
    elif lib == "cuda":
        initialize = _cuda_initialize
        random_kfloat = xoroshiro128p_uniform_float32
    else:
        raise NotImplementedError(f"Array library '{lib}' not currently handled")