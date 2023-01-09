# vect.py
"""Define Particles class and related support functions"""

from dataclasses import dataclass
import numpy as np

import torch

FLOAT_ARRAY_TYPE = np.ndarray
FLOAT_ARRAY_DTYPE = np.float32

INT_ARRAY_TYPE = np.ndarray
INT_ARRAY_DTYPE = np.int32
ZEROS_FN = np.zeros


class Particles:
    NUM_FLOAT_PARAMS = 11
    ENERGY, X, Y, Z, U, V, W, WT, DNEAR, TSTEP, USTEP = range(NUM_FLOAT_PARAMS)
    # For iparticles array (particle int parameters)
    NUM_INT_PARAMS = 3
    MEDIUM, REGION, STATUS = range(NUM_INT_PARAMS)
    int_keys = ("MEDIUM", "REGION", "STATUS")
    float_keys = "ENERGY X Y Z U V W WT DNEAR TSTEP USTEP".split()

    def __init__ (self, f_arr=None, i_arr=None, num_particles=None, array_library=None, **kwargs):
        if f_arr is not None:
            self.f_arr = f_arr
            self.i_arr = i_arr  # XXX should check is also not None
        else:
            raise NotImplementedError("Need to set arrays yourself, for now")
            # kwargs specify PARAMS to initialize to values other than zero
            # self.f_arr = FLOAT_ARRAY_TYPE.zeros(
            #     (self.NUM_FLOAT_PARAMS, num_particles) , dtype=FLOAT_ARRAY_DTYPE
            # )
            # self.i_arr = INT_ARRAY_TYPE.zeros(
            #     ( self.NUM_INT_PARAMS, num_particles), dtype=INT_ARRAY_DTYPE
            # )
            # self.set(**kwargs)
        self.array_libary = array_library

    def __len__(self):
        return len(self.f_arr[0])

    def any_alive(self):
        if self.array_libary == "pytorch":
            anyfunc = torch.any
        elif self.array_libary== "numpy":
            anyfunc = np.any
        else:
            raise NotImplementedError("Unknown array library")

        return anyfunc(self.energy > 0)

    def log_energy(self):
        return torch.log(self.energy)

    def keep(self, mask):
        """Return Particles with specified items according to mask
        """
        # if mask is all True, then just return self - no change.
        if all(mask):
            return self
        return Particles(self.f_arr[:, mask], self.i_arr[:, mask])

    def set(self, **kwargs):
        for key, val in kwargs.items():
            key = key.upper()
            index = getattr(self, key)
            if key in self.int_keys:
                self.i_arr[index] = val
            else:
                self.f_arr[index] = val

        # Return new Particles for GPU - can't modify in-place?
        return self
    def to_gpu(self):
        self.f_arr = self.f_arr.to("cuda")
        self.i_arr = self.i_arr.to("cuda")

    # Float array properties ---------------------
    @property
    def energy(self):
        return self.f_arr[self.ENERGY]

    @property
    def w(self):
        return self.f_arr[self.W]

    @property
    def z(self):
        return self.f_arr[self.Z]

    # Integer array properties ---------------------
    @property
    def medium(self):
        return self.i_arr[self.MEDIUM]

    @property
    def region(self):
        return self.i_arr[self.REGION]

    @property
    def status(self):
        return self.i_arr[self.STATUS]


if __name__ == "__main__":
    # Below is for testing code during dev...

    P = Particles  # short-form for accessing indices
    # Set up toy slab geometry like tutor examples
    num_particles = 10
    f_arr = ZEROS_FN(
                (P.NUM_FLOAT_PARAMS, num_particles) , dtype=FLOAT_ARRAY_DTYPE
            )
    i_arr = ZEROS_FN(
                (P.NUM_INT_PARAMS, num_particles), dtype=INT_ARRAY_DTYPE
            )


    # Start in region 0, medium 0, status 0 so leave those alone
    f_arr[P.ENERGY, :] = 100.0
    f_arr[P.W, :] = 1.0

    particles = Particles(f_arr, i_arr)
    print(particles.f_arr)
    print(particles.i_arr)
    z_bound = 100 # cm