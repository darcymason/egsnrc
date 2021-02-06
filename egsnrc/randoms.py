from egsnrc.egsfortran import ranlux, randomm
from egsnrc.commons import rng_seed, rng_array


def randomset():
    global rng_seed

    if rng_seed > 24:
        ranlux(rng_array)
        randomm.rng_seed = 1

    random_num = rng_array[rng_seed-1]
    randomm.rng_seed += 1

    return random_num


def random_iter():
    while True:
        yield from rng_array
        ranlux(rng_array)
