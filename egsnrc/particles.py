# particles.py
"""Define Particles class and related support functions"""

from dataclasses import dataclass
from collections import namedtuple
import numpy as np


particle_iattrs = "status region".split()
particle_fattrs = "energy x y z u v w".split()


STATUS, REGION = np.arange(2, dtype=np.int32)
ENERGY, Z = np.arange(2, dtype=np.int32)

# Constants to use for `status`
STEPPING, COMPTON, PHOTO = np.arange(3, dtype=np.int32)
# Negative ones for reason particle complete
GEOMETRY_DISCARD, PCUT_DISCARD, PHOTO_DISCARD, NO_DISCARD = np.arange(-3, 1, dtype=np.int32)

# PosDir = namedtuple("PosDir", ("x", "y", "z", "u", "v", "w"))
Particle = namedtuple(
    "Particle", (
        "status", "region", "energy", "x", "y", "z", "u", "v", "w"
    )
)

# type below not currently used
# ParticleType = NamedTuple((int32, int32, float32, float32), Particle)
