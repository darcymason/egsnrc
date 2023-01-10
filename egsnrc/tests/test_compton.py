import pytest
from egsnrc.expts.py.pycompton import py_compton  # XXX later need to switch to live

from egsnrc import random

def test_expt_pycompton_1MeV():
    """Test known outputs from modified 1 MeV photon tutor7 (mortran)"""
    # NOTE: have to not do random azimuth angle in the test compton code to
    #    avoid consuming randoms not captured in tutor7. Pass argument to func.
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

def test_expt_pycompton_10MeV():
    """Test known outputs from modified 10 MeV photon tutor7 (mortran)"""
    # Set up list of known randoms from original tutor7 run
    random.set_array_library("sequence")
    rand_sequence = ([
        0.42048686742782593, 0.49585634469985962, 0.57738614082336426,
        7.6826930046081543E-002, 0.40864366292953491, 0.97492825984954834,
        0.88091498613357544, 0.96854400634765625, 0.17157137393951416,
        0.52135562896728516, 0.11907631158828735, 0.48996019363403320,
        0.94912832975387573, 0.99772697687149048,
        0.41744160652160645, 8.0652594566345215E-002,
        0.10784226655960083, 4.6389222145080566E-002,
        2.8621673583984375E-002, 8.8309347629547119E-002, 0.19784265756607056,
        0.81320518255233765, 0.62365549802780151,
        0.32891792058944702, 0.41770315170288086,
    ])
    rng = random.initialize(rand_sequence)

    # Energies from tutor7 run with tutor7_10MeV_Si.egsinp file
    # Had multiple interactions before exiting so this tests various energies
    energies = [
        10.0, 1.5544328853911626, 10.0, 0.38670521268561409, 0.25105326292747948,
        10.0, 7.3636066316550799E-002, 0.34517802901394307
    ]
    result_e_costhe = [  # NOTE this is "costhe before changes" i.e. before UPHI etc.
        (1.5544328853911626, 0.72236330760030842),
        (1.5302843861169189, 0.99481240523080916),
        (0.38670521268561409, -0.27031741977450796),
        (0.25105326292747948, 0.28599679089128460),
        (0.14004668573376408, -0.61335532134259729),
        (0.34517802901394307, -0.42929260065073582),
        (7.0558756714006632E-002, 0.69734370445616700),
        (0.21206500243906989, 7.0758854413525252E-002),
    ]

    for e_in, e_costhe in zip(energies, result_e_costhe):
        result = py_compton(rng, e_in, calc_azimuth=False)
        energy, _, costhe = result[:3]
        assert (energy, costhe) == pytest.approx(e_costhe)
