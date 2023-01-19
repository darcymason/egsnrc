# example1.py
"""Define a slab geometry with photon tracking only
"""
import numba as nb
from numba import cuda
import numpy as np

from egsnrc.photon import COMPTON, PHOTO_DISCARD
from egsnrc.particles import Particle

NO_DISCARD, DISCARD = np.arange(2, dtype=np.int32)

have_gpu = cuda.is_available()

SCORE_nCOMPTON, SCORE_nPHOTO, SCORE_eCOMPTON, SCORE_ePHOTO, \
    SCORE_nLost, SCORE_eLost = np.arange(6, dtype=np.int32)

# Example of scoring function (ausgab) - replace with user one

def scoring(gid, p1, p2, out):
    # Note, if don't track by gid, then need "atomic" operations to handle multi-thread
    if p1.status == COMPTON:
        out[gid, p2.region, SCORE_nCOMPTON] += 1.0
        out[gid, p2.region, SCORE_eCOMPTON] += p1.energy - p2.energy
    elif p1.status == PHOTO_DISCARD:
        out[gid, p2.region, SCORE_nPHOTO] += 1.0
        out[gid, p2.region, SCORE_ePHOTO] += p1.energy - p2.energy
    elif p2.region == NUM_REGIONS:  # lost to geometry (0-based region index)
        out[gid, p2.region, SCORE_nLost] += 1.0
        out[gid, p2.region, SCORE_eLost] += p1.energy

if have_gpu:
    scoring = cuda.jit(scoring, device=True)  # "Device function"


boundaries = np.array([-99.0, 0.0, 2.0, 4.0], dtype=np.float32)

def howfar(p, ustep):  # -> step, region, discard_flag (>0 to discard)
    """Given particle and proposed step distance, return actual step and region

    If crossing a region boundary, step is distance to the new region,
    and region returned is new region.
    """
    # Slab geometry with "z" going to right:
    # Region 0       |     Region 1    |    Region 2      |    Region 3
    #   vacuum       |                 |                  |    vacuum
    #                |     Medium0     |     Medium1      |
    #       photon-> |                 |                  |
    #                |                 |                  |
    #               z=0              2 cm               4 cm


    if p.region not in (1, 2):  # outside the geometry
        return ustep, p.region, DISCARD

    # Else in a slab
    # Deal with trivial case first, simplifies the other cases into one
    if p.w == 0.0:  # parallel to slabs, can't cross any
        return ustep, p.region, NO_DISCARD

    boundary_offset = 1 if p.w > 0.0 else 0
    region_change = 1 if p.w > 0.0 else -1

    tval = (boundaries[p.region + boundary_offset] - p.z) / p.w
    if tval > ustep:  # step doesn't reach boundary
        return ustep, p.region, NO_DISCARD

    # Return step to just reach the boundary and set new region
    return tval, p.region + region_change, NO_DISCARD


if have_gpu:
    howfar = cuda.jit(howfar, device=True)


def test_gpu():
    iparticles = np.zeros((5, 2), dtype=np.int32)
    # set regions
    iparticles[:, 1] = [1, 1, 1, 2, 2]
    fparticles = np.zeros((5, 7), dtype=np.float32)
    fparticles[:, 3] = [0.5, 1.5, 1.75, 2.2, 3.8]  # z values
    fparticles[:, -1] = [1, -1, 0, -0.8, -0.5] # w values
    usteps = np.array([5, 3, 2, 1, 8], dtype=np.float32)
    fout = np.zeros(len(fparticles), dtype=np.float32)
    iout = np.zeros((len(iparticles), 2), dtype=np.int32)
    test_kernel.forall(len(iparticles))(
        usteps, iparticles, fparticles, iout, fout
    )
    cuda.synchronize()
    print(iout)
    print(fout)
    # [[2 0]
    #  [0 0]
    #  [1 0]
    #  [1 0]
    #  [1 0]]
    # [1.5        1.5        2.         0.25000006 3.6       ]

# Dummy kernel just for testing howfar
@cuda.jit
def test_kernel(usteps, iparticles, fparticles, iout, fout):
    i = cuda.grid(1)
    if i > len(iparticles):
        return
    ip = iparticles[i]
    fp = fparticles[i]
    p = Particle(ip[0], ip[1], fp[0], fp[1], fp[2], fp[3], fp[4], fp[5], fp[6])
    step, region, discard = howfar(p, usteps[i])
    fout[i] = step
    iout[i, 0] = region
    iout[i, 1] = discard


if __name__ == "__main__":
    pos = (0, 0, 1.7)
    dir = (0.0, 0.0, -1)
    p = Particle(0, 1, 1.0, *pos, *dir)
    print(f"{p=}")
    print(f"{howfar(p, 1.5)=}")
    print(f"{howfar(p, 5.5)=}")
    # test_gpu()
