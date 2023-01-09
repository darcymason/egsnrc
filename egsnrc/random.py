
import numpy as np
try:
    import torch
except ImportError:
    pass

# Numpy random functions
def _np_initialize(seed):
    return np.random.default_rng(seed)


def _np_floats_0_1(rng, num, device=None):
    return rng, rng.random(num)


# PyTorch random funtions ---------
def _torch_initialize(seed):
    return torch.random.manual_seed(seed)


def _torch_floats_0_1(key, num, device=None):
    return key, torch.rand(num, device=device)


# Configuration to select between libraries
def set_array_library(lib: str):
    global floats_0_1, initialize

    if lib == "numpy":
        floats_0_1 = _np_floats_0_1  # Later switch depending on using gpu arrays
        initialize = _np_initialize
    elif lib == "pytorch":
        floats_0_1 = _torch_floats_0_1
        initialize = _torch_initialize
    else:
        raise NotImplementedError("Array library not currently handled")