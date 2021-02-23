import pytest
from egsnrc.util import float_from_fort_hex

def test_py_hex_from_fort_hex():
    vals = [
        ("            0", 0.0),
        ("3FF3C083126E978D", 1.2345),
        ("BFEFAE147AE147AE",  -0.99),
        ("3FD55553EF6B5D46",  0.333333),
        ("BFD55553EF6B5D46", -0.333333),
        ("8000000000000065", -5e-322),
        # couple of 'denormal' numbers
        ("           65",  5e-322),
        ("            1",  5e-324),
    ]

    for fort_hex, expected in vals:
        assert expected == float_from_fort_hex(fort_hex)
