
import itertools
import numpy as np
try:
    import torch
except ImportError:
    pass

# Numpy random functions
def _np_initialize(seed):
    return np.random.default_rng(seed)


def _np_floats_0_1(rng, num, device=None):
    rands = rng.random(num)
    # print("Random floats:", rands)
    return rng, rands


# PyTorch random funtions ---------
def _torch_initialize(seed):
    return torch.random.manual_seed(seed)


def _torch_floats_0_1(key, num, device=None):
    return key, torch.rand(num, device=device)

# Set up "random" from known sequence -------------
#    (e.g. to compare with mortran EGSnrc)
def _seq_initialize(known_list, vect=False, list_type=np.array):
    return SeqGen(known_list, vect, list_type)

def _seq_floats_0_1(rng, num, device=None):
    return (rng, *rng.random(num))

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
    global floats_0_1, initialize

    if lib == "numpy":
        floats_0_1 = _np_floats_0_1  # Later switch depending on using gpu arrays
        initialize = _np_initialize
    elif lib == "pytorch":
        floats_0_1 = _torch_floats_0_1
        initialize = _torch_initialize
    elif lib == "sequence":
        floats_0_1 = _seq_floats_0_1
        initialize = _seq_initialize
    else:
        raise NotImplementedError("Array library not currently handled")