import pytest
from egsnrc.expts.py.pycompton import py_compton  # XXX later need to switch to live

from egsnrc import random

def test_expt_pycompton():
    """Test known outputs from modified 1 MeV photon tutor7 (mortran)"""
    # NOTE: have to comment out the random azimuth angle in the test compton code
    random.set_array_library("sequence")
    rand_sequence = ([
        0.22987532615661621, 0.78682494163513184,
        0.76323693990707397, 0.19077306985855103,
        0.55831658840179443,  4.9074590206146240E-002,
        0.55832588672637939, 0.29218196868896484,
        0.65860986709594727,  4.0171027183532715E-002,
        1.7512619495391846E-002,  0.44601458311080933,
        0.91792672872543335, 0.54524618387222290,
        0.28818345069885254,  3.5068333148956299E-002,
        0.13349604606628418, 0.85515218973159790,
        0.54984831809997559,  5.7695209980010986E-002,
        0.49270433187484741, 0.94706720113754272,
        0.84880810976028442, 0.67912405729293823,
        0.30815130472183228, 0.37652158737182617,
        0.96473675966262817, 0.67000657320022583,
        0.96202033758163452, 0.26576608419418335,
        0.62973368167877197, 0.13346099853515625,
    ])
    energies = [1.0] * 10
    result_e_costhe = [  # NOTE this is "costhe before changes" i.e. before UPHI etc.
        (0.81141922919270071, 0.88123946893996008),
        (0.64820104040175874, 0.72266489438445403),
        (0.64820844647776044, 0.72267390145936683),
        (0.72808421050232852, 0.80915849456905908),
        (0.21745297830124174, -0.83892959634582454),
        (0.43304114464512694, 0.33097491664454248),
        (0.64145609910772927, 0.71437552426368434),
        (0.96974936034766024, 0.98405975176254867),
        (0.70508445157959732, 0.78626455389275030),
    ]

    rng = random.initialize(rand_sequence)
    for e_in, e_costhe in zip(energies, result_e_costhe):
        result = py_compton(rng, e_in, calc_azimuth=False)
        energy, _, costhe = result[:3]
        assert (energy, costhe) == pytest.approx(e_costhe)
