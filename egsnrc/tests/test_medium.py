import pytest
from egsnrc.hatch import DATA_DIR, get_xsection_table
from egsnrc.media import Medium, Element, Interaction


class TestSigma:
    def test_one_element_photo(self):
        """Simple Medium with one element calculates correct sigma"""
        PHOTO = Interaction.PHOTOELECTRIC
        photo_data = get_xsection_table(DATA_DIR / "xcom_photo.data")
        element_ta = Element(z=73, pz=1, wa=0)
        medium_ta = Medium("Ta", [element_ta], rho=1, ap=0.01, up=50.0)
        medium_ta.calc_sigma(PHOTO, photo_data)
        print(medium.sigma[PHOTO])
        assert False