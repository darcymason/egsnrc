# regions.py
"""Define Region class and related support functions"""

from dataclasses import dataclass
from collections import namedtuple
import numpy as np
from egsnrc.config import device_jit
from egsnrc.media import _Medium
from numba import cuda
import numpy as np


MEDIUM, IRAYL, IPHOTONUCR = np.arange(3, dtype=np.int32)
PCUT, RHO = np.arange(2, dtype=np.float32)


_Region = namedtuple("Region", ("number medium pcut rho"))  # irayl iphotonucr

# @dataclass
# class Region:
#     region_number: int
#     medium: Medium
#     pcut: float
#     rho: float

# # Internal function to convert classes to Region namedtuples to pass to kernels
# def _regions(regions, media):
#     return tuple(
#         _Region(r.region_number, r.medium.)
#     )


# @device_jit
# def set_region(i, iregions, fregions):
#     """Return a Region namedtuple for the given array-based info

#     Used in e.g. photon kernel
#     """
#     oi = iregions[i]
#     of = fregions[i]

#     return Region(oi[MEDIUM], of[PCUT], of[RHO])  # oi[IRAYL], oi[IPHOTONUCR],

def Region(number, medium, pcut=0.001, rho=None):
    if not isinstance(medium, _Medium):
        raise TypeError("medium must be a Medium named tuple")
    if rho is None:
        rho = medium.rho
    return _Region(np.int32(number), medium, np.float32(pcut), np.float32(rho))
