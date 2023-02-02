# particles.py
"""Define Particle object and related support functions"""

from collections import namedtuple
import numpy as np
from egsnrc.config import device_jit
from egsnrc.params import EPSGMFP, VACDST


# Constants to use for `status`
STEPPING, COMPTON, PHOTO = np.arange(3, dtype=np.int32)
# Negative ones for reason particle complete
INTERACTION_READY, GEOMETRY_DISCARD, PCUT_DISCARD, USER_PHOTON_DISCARD, \
    PHOTO_DISCARD, NO_DISCARD = np.arange(-5, 1, dtype=np.int32)


# Particle attributes
particle_iattrs = "status region".split()
particle_fattrs = "energy x y z u v w".split()

STATUS, REGION = np.arange(2, dtype=np.int32)
ENERGY, X, Y, Z, U, V, W = np.arange(7, dtype=np.int32)

Particle = namedtuple("Particle", "status region energy x y z u v w")

# Functions to manage particles here, because with namedtuples it is a little
#  trickier to remember order, and this allows all changes in one place
#  if attributes are later added to the namedtuple
@device_jit
def set_particle(i, regions, iparticles, fparticles):
    """Return a Particle for the given array-based info

    Used in e.g. photon kernel
    """
    oi = iparticles[i]
    of = fparticles[i]
    return Particle(
        oi[STATUS], regions[oi[REGION]],
        of[ENERGY], of[X], of[Y], of[Z], of[U], of[V], of[W])


@device_jit
def replace_e_uvw(p, energy, u, v, w):
    """Return a Particle like `p` but with uvw replaced"""
    return Particle(p.status, p.region, energy, p.x, p.y, p.z, u, v, w)

@device_jit
def replace_region_xyz(p, region, x, y, z):
    """Return a Particle like `p` but with x,y,z replaced"""
    return Particle(p.status, region, p.energy, x, y, z, p.u, p.v, p.w)


class PhotonSource:
    """Convenience class to define a photon source

    Currently can specify a min and max energy to sample randomly from,
    but only a fixed position and direction

    Internally, supplies the types needed to pass to the GPU kernel
    """
    # Could easily add more to this class,
    # e.g. uniform direction with a `direction=None` default
    #  and line sources, etc. with multi-dimensional `position`
    # Could pass callbacks for any of these for user-specified spectra, etc.
    def __init__(self, energy, region, position, direction, store_particles=False):
        if isinstance(energy, (tuple, list)) and len(energy) != 2:
            raise ValueError("energy must be a single value or two-tuple")
        self.energy = energy
        self.region = region
        self.position = position
        self.direction = direction
        self.total_energy = 0
        self.store_particles = store_particles
        if self.store_particles:
            self.iparticles = np.empty((0, len(particle_iattrs)), dtype=np.int32)
            self.fparticles = np.empty((0, len(particle_fattrs)), dtype=np.float32)

    def generate(self, num_particles):
        # Start with zeros, then fill in
        iparticles = np.zeros(
            (num_particles, len(particle_iattrs)), dtype=np.int32
        )
        fparticles = np.zeros(
            (num_particles, len(particle_fattrs)), dtype=np.float32
        )

        energy = self.energy
        # Set requested energies
        if isinstance(energy, (tuple, list)):
            # XXX need to make into proper GPU kernel for generating random energies
            np.random.seed(42) # XXX
            energies = np.random.random(num_particles)*(energy[1] - energy[0]) + energy[0]
            energies = energies.astype(np.float32)
            fparticles[:, ENERGY] = energies
        else:
            fparticles[:, ENERGY] = energy

        self.total_energy += sum(fparticles[:, ENERGY])
        # Set region and position, direction
        iparticles[:, REGION] = self.region.number
        fparticles[:, X: X + 3] = self.position
        fparticles[:, U: U + 3] = self.direction

        if self.store_particles:
            self.fparticles = np.vstack((self.fparticles, fparticles))
            self.iparticles = np.vstack((self.iparticles, iparticles))
        return iparticles, fparticles


if __name__ == "__main__":
    source = PhotonSource(40000, (0.511, 1.511), 1, (0, 0, 0.5), (0, 0, 1))
    print(source.iparticles[:10])
    print(source.fparticles[:10])
    min_e = np.min(source.fparticles[:, ENERGY])
    max_e = np.max(source.fparticles[:, ENERGY])
    print(f"min, max e: {min_e}, {max_e}")
