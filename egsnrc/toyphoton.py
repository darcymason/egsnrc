import numpy as np
import torch.random
import torch
# =====================================
# RANDOM functions
import numpy.random as np_random
def _np_initialize(seed):
    np_random.seed(seed)
    return seed

def _np_floats_0_1(key, num):
    return key, np_random.random(num)

def _torch_initialize(seed):
    return torch.random.manual_seed(seed)

def _torch_floats_0_1(key, num):
    return key, torch.rand(num)


# ====================================================
# PARTICLES defn

class Particles:
    NUM_FLOAT_PARAMS = 11
    ENERGY, X, Y, Z, U, V, W, WT, DNEAR, TSTEP, USTEP = range(NUM_FLOAT_PARAMS)
    # For iparticles array (particle int parameters)
    NUM_INT_PARAMS = 3
    MEDIUM, REGION, STATUS = range(NUM_INT_PARAMS)
    int_keys = ("MEDIUM", "REGION", "STATUS")
    float_keys = "ENERGY X Y Z U V W WT DNEAR TSTEP USTEP".split()

    def __init__ (self, f_arr=None, i_arr=None, num_particles=None, **kwargs):
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

    def __len__(self):
        return len(self.f_arr[0])

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

# Set up which path, numpy or jax, etc.
def set_array_library(library="numpy"):
    global random_floats_0_1, random_initialize
    global FLOAT_ARRAY_TYPE, FLOAT_ARRAY_DTYPE
    global INT_ARRAY_TYPE, INT_ARRAY_DTYPE
    global ZEROS_FN

    if library == "numpy":
        random_floats_0_1 = _np_floats_0_1
        random_initialize = _np_initialize
        FLOAT_ARRAY_TYPE = np.ndarray
        FLOAT_ARRAY_DTYPE = np.float32

        INT_ARRAY_TYPE = np.ndarray
        INT_ARRAY_DTYPE = np.int32
        ZEROS_FN = np.zeros
    elif library == "pytorch":
        random_floats_0_1 = _torch_floats_0_1
        random_initialize = _torch_initialize
        FLOAT_ARRAY_TYPE = torch.Tensor
        FLOAT_ARRAY_DTYPE = torch.float32

        INT_ARRAY_TYPE = torch.Tensor
        INT_ARRAY_DTYPE = torch.int32
        ZEROS_FN = torch.zeros
    else:
        raise ValueError("Unknown library argument")
# ======================================
# PHOTON "transport"
def toy(particles, key):
    key, ran_floats = random_floats_0_1(key, len(particles.energy))
    energy = particles.energy - 20.0 * ran_floats
    pos_e_mask = energy > 0
    particles = particles.keep(pos_e_mask)
    particles = particles.set(energy=energy[pos_e_mask])

    return particles, key


# =========--------================-----------================---------========
# MAIN CODE

# Set up which path, numpy or jax, etc.
def set_array_library(library="numpy"):
    global random_floats_0_1, random_initialize
    global FLOAT_ARRAY_TYPE, FLOAT_ARRAY_DTYPE
    global INT_ARRAY_TYPE, INT_ARRAY_DTYPE
    global ZEROS_FN

    if library == "numpy":
        random_floats_0_1 = _np_floats_0_1
        random_initialize = _np_initialize
        FLOAT_ARRAY_TYPE = np.ndarray
        FLOAT_ARRAY_DTYPE = np.float32

        INT_ARRAY_TYPE = np.ndarray
        INT_ARRAY_DTYPE = np.int32
        ZEROS_FN = np.zeros
    elif library == "pytorch":
        random_floats_0_1 = _torch_floats_0_1
        random_initialize = _torch_initialize
        FLOAT_ARRAY_TYPE = torch.Tensor
        FLOAT_ARRAY_DTYPE = torch.float32

        INT_ARRAY_TYPE = torch.Tensor
        INT_ARRAY_DTYPE = torch.int32
        ZEROS_FN = torch.zeros
    else:
        raise ValueError("Unknown library argument")

def main(num_particles, array_library, gpu_requested=False):
    have_cuda = torch.cuda.is_available()
    if gpu_requested and (not have_cuda or array_library != "pytorch"):
        raise ValueError("GPU not available")

    set_array_library(array_library)

    P = Particles  # short-form for accessing indices
    # Set up toy slab geometry like tutor examples

    f_arr = ZEROS_FN(
                (P.NUM_FLOAT_PARAMS, num_particles) , dtype=FLOAT_ARRAY_DTYPE
            )
    i_arr = ZEROS_FN(
                (P.NUM_INT_PARAMS, num_particles), dtype=INT_ARRAY_DTYPE
            )

    key = random_initialize(42)

    # Start in region 0, medium 0, status 0 so leave those alone
    f_arr[P.ENERGY, :] = 100.0
    f_arr[P.W, :] = 1.0

    particles = Particles(f_arr, i_arr)
    if have_cuda:
        particles.to_gpu()

    # print(particles.f_arr)
    # print(particles.i_arr)
    z_bound = 100 # cm

    while len(particles):
        # print(f"{len(particles)=}")
        particles, key = toy(particles, key)

import timeit


num_particles = 10_000
lib = "numpy"
gpu = False  # True
lib = "pytorch"

gpu_msg = f" with GPU" if gpu else "** NO GPU **"
print(f"Starting run(s) {gpu_msg}...")
seconds = timeit.timeit("main(num_particles, lib, gpu)", globals=globals(), number=10)
print(f"Completed run in {seconds} seconds")
# main(10_000, 'numpy')


def setup(num_particles, array_library, gpu_requested=False):
    have_cuda = torch.cuda.is_available()
    device = "cuda" if (have_cuda and gpu_requested) else "cpu"
    if gpu_requested and (not have_cuda or array_library != "pytorch"):
        raise ValueError("GPU not available")

    set_array_library(array_library)

    P = Particles  # short-form for accessing indices
    # Set up toy slab geometry like tutor examples

    f_arr = ZEROS_FN(
                (P.NUM_FLOAT_PARAMS, num_particles) , dtype=FLOAT_ARRAY_DTYPE
            )
    i_arr = ZEROS_FN(
                (P.NUM_INT_PARAMS, num_particles), dtype=INT_ARRAY_DTYPE
            )

    key = random_initialize(42)

    # Start in region 0, medium 0, status 0 so leave those alone
    f_arr[P.ENERGY, :] = 100.0
    f_arr[P.W, :] = 1.0

    particles = Particles(f_arr, i_arr)
    if have_cuda and gpu_requested:
        particles.to_gpu()

    # print(particles.f_arr)
    # print(particles.i_arr)
    z_bound = 100 # cm