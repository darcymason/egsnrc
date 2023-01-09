import numpy as np
import pytest
from egsnrc.hatch import DATA_DIR, get_xsection_table
from egsnrc.media import Medium, Element, Interaction
from math import exp, log



def test_sigmas_photo(self):
    """Simple Medium with one element calculates correct sigmas"""
    PHOTO = Interaction.PHOTOELECTRIC
    photo_data = get_xsection_table(DATA_DIR / "xcom_photo.data")
    element = Element(z=73, pz=1)

    # Set medium lower E to same as photo_data and num steps to double the input data
    first_energy = exp(photo_data[73][0][0])
    first_sigma = exp(photo_data[73][1][0])
    medium = Medium(
        "Ta", [element], rho=1, ap=first_energy, up=50.0, mge=2 * len(photo_data)
    )
    # medium._calc_sigmas(PHOTO, photo_data) now done in Medium.__post_init__
    sigmas = medium.sigmas[PHOTO]
    # Trivial test - first sigma as expected
    assert sigmas[0] == first_sigma

    # Now interpolate more using high-res num steps and different AP than table
    element = Element(z=6, pz=1)
    medium = Medium("C", [element], rho=1, ap=1.0, up=25.0, mge=5_000)
    # medium._calc_sigmas(PHOTO, photo_data) now done in Medium.__post_init__

    # photo table for Carbon excerpted here around point of interest
    # 2.89037 -14.44193   2.99573 -14.55324   3.09104 -14.65345   3.17805 -14.74463
    test_E = exp(2.99573)
    expect_sigma = exp(-14.55324)

    assert medium.calc_sigma(PHOTO, test_E) == pytest.approx(expect_sigma)
