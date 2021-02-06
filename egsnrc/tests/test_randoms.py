import pytest

pytest.importorskip("egsnrc.egsfortran")  # from numpy.f2py, used while in transition
from egsnrc import egsfortran
from egsnrc.randoms import randomset, random_iter, rng_array


class TestRandom:
    def test_egsfor_ranlux_randomset(self):
        """Test RANLUX random number generator"""
        # Based on ranlux_test.mortran in EGSnrc
        # Assumes default egsfortran random macros
        egsfortran.init_ranlux(1, 1)
        egsfortran.ranlux(rng_array)
        egsfortran.randomm.rng_seed = 1
        first10_sum = sum(randomset() for x in range(10))
        assert 6.041212022304535 == first10_sum
        first10k_sum = first10_sum + sum(randomset() for x in range(9990))
        assert 5037.366532325745 == first10k_sum
        # Add back first million if desired - leave out to save a few seconds
        # first_million_sum = first10k_sum + sum(randomset() for x in range(990_000))
        # assert 500181.8234628493 == first_million_sum

    def test_egsfor_ranlux_iter(self):
        """Test RANLUX random number generator as an iterator"""
        # Based on ranlux_test.mortran in EGSnrc
        # Assumes default egsfortran random macros
        egsfortran.init_ranlux(1, 1)
        egsfortran.ranlux(rng_array)
        rnd_iter = random_iter()
        first10_sum = sum(next(rnd_iter) for i in range(10))
        assert 6.041212022304535 == first10_sum
        first10k_sum = first10_sum + sum(next(rnd_iter) for x in range(9990))
        assert 5037.366532325745 == first10k_sum
        first_million_sum = first10k_sum + sum(next(rnd_iter) for x in range(990_000))
        assert 500181.8234628493 == first_million_sum

