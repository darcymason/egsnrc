
import itertools
import numpy as np
import numba

from numba.cuda.random import create_xoroshiro128p_states
from numba.cuda.random import xoroshiro128p_uniform_float32


def _np_initialize(seed, num_particles=None):
    return np.random.default_rng(seed)

def _np_float32(rng, i):
    return rng.random()

# Set up "random" from known sequence -------------
#    (e.g. to compare with mortran EGSnrc)
def _seq_initialize(known_list, vect=False, list_type=np.array):
    return SeqGen(known_list, vect, list_type)

def _seq_float32(rng, num, device=None):
    return (rng, *rng.random(num))

def _cuda_initialize(seed, num_particles):
    return numba.cuda.to_device(create_xoroshiro128p_states(num_particles, seed))

class SeqGen:
    def __init__(self, known_list, vect=False, list_type=np.array):
        # list_type could be e.g. torch.Tensor
        if not vect:
            self.known_list = list_type(
                itertools.chain.from_iterable(
                    itertools.chain.from_iterable(known_list)
                )
            )
        else:
            # Transpose the list so each "index" in the vector gets its
            #    proper original random number sequence
            self.known_list = list_type(list(
                itertools.chain.from_iterable(
                    itertools.chain.from_iterable(
                        filter(None, sublist)
                        for sublist in itertools.zip_longest(*known_list)
                    )
                )
            )
            )
        self.pos = 0
    def random(self, num):
        shape = False
        if isinstance(num, tuple):
            shape = num
            num = shape[0] * shape[1]
        self.pos += num
        rands = self.known_list[self.pos - num: self.pos]
        if shape:
            rands = rands.reshape((shape[1], shape[0]))
            # print(rands)
            return (*(rands[:,i] for i in range(shape[0])),)
        # rands = rands.transpose()
        return rands


# Configuration to select between libraries
def set_array_library(lib: str):
    global random_float32
    global initialize

    if lib == "numpy":
        random_float32 = _np_float32  # Later switch depending on using gpu arrays
        initialize = _np_initialize
    elif lib == "sequence":
        random_float32 = _seq_float32
        initialize = _seq_initialize
    elif lib == "cuda":
        initialize = _cuda_initialize
        random_float32 = xoroshiro128p_uniform_float32
    else:
        raise NotImplementedError(f"Array library '{lib}' not currently handled")