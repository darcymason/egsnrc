# example1.py
"""Define a slab geometry with photon tracking only
"""
# import os
# os.environ['NUMBA_ENABLE_CUDASIM'] = '1'  # XXX temp for testing GPU running on CPU

import numba as nb
from numba import cuda
import numpy as np

from egsnrc.photon import COMPTON, PHOTO
from egsnrc.particles import (
    Particle, particle_iattrs, particle_fattrs, ENERGY, W, REGION, PhotonSource
)
from egsnrc.config import device_jit
from egsnrc.media import Medium, Vacuum
from egsnrc.regions import Region
from egsnrc.shower import shower


NO_DISCARD, DISCARD = np.arange(2, dtype=np.int32)


SCORE_nCOMPTON, SCORE_nPHOTO, SCORE_nLOST = np.arange(3, dtype=np.int32)
SCORE_eCOMPTON, SCORE_ePHOTO, SCORE_eLOST = np.arange(3, dtype=np.int32)


# Example of scoring function (ausgab)
@device_jit
def ausgab(gid, status, p1, p2, iscore, fscore):
    # Note, if don't track by gid, then need "atomic" operations to handle multi-thread
    region_number = p2.region.number
    if status == COMPTON:
        iscore[gid, region_number, SCORE_nCOMPTON] += 1
        fscore[gid, region_number, SCORE_eCOMPTON] += p1.energy - p2.energy
    elif status == PHOTO:
        iscore[gid, region_number, SCORE_nPHOTO] += 1
        fscore[gid, region_number, SCORE_ePHOTO] += p1.energy - p2.energy
    elif region_number not in (1, 2):  # lost to geometry
        iscore[gid, region_number, SCORE_nLOST] += 1
        fscore[gid, region_number, SCORE_eLOST] += p1.energy


@device_jit
def howfar(p, regions, ustep):  # -> step, region, discard_flag (>0 to discard)
    """Given particle and proposed step distance, return actual step and region

    If crossing a region boundary, step is distance to the new region,
    and region returned is new region.
    """
    # Slab geometry with "z" going to right:
    # Region 0       |     Region 1    |    Region 2      |    Region 3
    #   vacuum       |                 |                  |    vacuum
    #                |     Medium1     |     Medium2      |
    #         photon |->               |                  |
    #                |                 |                  |
    #               z=0              2 cm               4 cm


    region = p.region
    region_num = region.number
    if region_num not in (1, 2):  # outside the geometry
        return ustep, region, DISCARD

    # Else in a slab
    # Deal with trivial case first, simplifies the other cases into one
    if p.w == 0.0:  # parallel to slabs, can't cross any
        return ustep, region, NO_DISCARD

    boundary_offset = 1 if p.w > 0.0 else 0
    region_change = 1 if p.w > 0.0 else -1

    tval = (boundaries[region_num + boundary_offset] - p.z) / p.w
    if tval > ustep:  # step doesn't reach boundary
        return ustep, region, NO_DISCARD

    # Return step to just reach the boundary and set new region
    return tval, regions[region_num + region_change], NO_DISCARD


usage = """
python example1.py [num_particles] [num_batches]

If optional num_particles is not specified, it defaults to 50 for testing.
If num_batches is not specified, it defaults to 1 on CPU or 2 on GPU (to separate compile time)
"""

if __name__ == "__main__":
    import sys
    num_particles = 12
    if len(sys.argv) > 1:
        try:
            num_particles = int(sys.argv[1])
        except ValueError:
            print("Optional num_particles command-line argument must be an integer")
            print(usage)
            sys.exit(-1)
    if len(sys.argv) > 2:
        try:
            num_batches = int(sys.argv[2])
        except ValueError:
            print("Optional num_batches command-line argument must be an integer")
            print(usage)


    # Set up the media and the regions
    Ta = Medium(1, "Ta")
    Si = Medium(2, "Si")
    media = [Ta, Si]

    regions = [
        Region(0, Vacuum, 0.001, 0),
        Region(1, Ta, 0.001, Ta.rho),
        Region(2, Si, 0.001, Si.rho),
        Region(3, Vacuum, 0.001, 0)
    ]


    # Boundaries of slabs (cm)
    # First one is dummy because of 0-based indexing and region 0 being outside geom
    boundaries = np.array([-1.0e8, 0.0, 2.0, 4.0], dtype=np.float32)

    # Set up particles with random range of energies
    # starting at (0, 0, 0) at edge of Region 1, traveling in z direction
    source = PhotonSource(
        energy=(0.511, 1.0), region=regions[1],
        position=(0, 0, 0), direction=(0, 0, 1),
        store_particles=True  # for debugging with small num_particles
    )
    # Initialize Scoring
    NUM_REGIONS = 4

    # Initialize Scoring arrays
    iscore = np.zeros(
        (num_particles, NUM_REGIONS,  3),  # (nCompt, nPhoto, nLost)
        dtype=np.int32
    )
    fscore = np.zeros(
        (num_particles, NUM_REGIONS,  3),  # (eCompt, ePhoto, eLost)
        dtype=np.float32
    )

    # Start the simulation
    shower(
        42, num_particles, source, regions, media, howfar, ausgab, iscore, fscore
    )


    if len(fscore) < 51:
        print("fparticle energies")
        print(source.fparticles[:, ENERGY])
        print("Scoring - Particles/Region")
        print("nCompt     nPhoto     eCompt     ePhoto     nLost      eLost")
        print(fscore)
    for r in range(1, 3):
        print("Energy Total by Region")
        print(f"Region {r}----")
        print("Compton: ", sum(fscore[:,r,SCORE_eCOMPTON]))
        print("Photo  : ", sum(fscore[:,r,SCORE_ePHOTO]))
    sum_compt = np.sum(fscore[:, :, SCORE_eCOMPTON])
    sum_photo = np.sum(fscore[:, :, SCORE_ePHOTO])
    sum_lost = np.sum(fscore[:, :, SCORE_eLOST])
    print(f"Sums:  Compt: {sum_compt}, Photo {sum_photo}, Lost {sum_lost}")
    print("Energy in :", source.total_energy)
    print("Energy out:", sum((sum_compt, sum_photo, sum_lost)))

