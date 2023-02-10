import numpy as np
import pytest
from egsnrc.config import KFLOAT
from egsnrc.hatch import DATA_DIR, get_xsection_table
from egsnrc.media import Medium, Interaction
from math import exp, log
from pathlib import Path

HERE = Path(__file__).resolve().parent


def test_sigmas_photo():
    """Simple Medium with one element calculates correct sigmas"""
    PHOTO = Interaction.PHOTOELECTRIC
    photo_data = get_xsection_table(DATA_DIR / "xcom_photo.data")

    first_energy = exp(photo_data[6][0][0])
    first_sigma = exp(photo_data[6][1][0])

    # Set medium lower E to same as photo_data and num steps to double the input data
    medium = Medium(1, "C", rho=1, ap=first_energy, up=50.0, mge=2 * len(photo_data))
    sigmas = medium.sigmas[PHOTO]
    # Trivial test - first sigma as expected
    assert sigmas[0] == pytest.approx(first_sigma)

    # Now interpolate more using high-res num steps and different AP than table
    medium = Medium(2, "C", ap=1.0, up=25.0, mge=5_000)

    # photo table for Carbon excerpted here around point of interest
    # 2.89037 -14.44193   2.99573 -14.55324   3.09104 -14.65345   3.17805 -14.74463
    test_E = exp(2.99573)
    expect_sigma = exp(-14.55324)

    assert medium.calc_sigma(PHOTO, test_E) == pytest.approx(expect_sigma)

def test_gmfp_and_gbr_vs_mortran():
    """Calculated gmfp's and branching ratios match known values from mortran"""
    # Si data produced by modified tutor7.mortran in
    # egsnrc/egsnrc/egs_home/samp_mfp_br
    # PARAMETERS $MXGE=20, which was output to gen_sigmas_sample_Si.xsections,
    # by modified subroutine `egs_init_user_photon` of egsnrc.mortran
    # outputs modified to give for i=1 to nge:
    # GMFP0   GMFP1   GBR10     GBR11    GBR20    GBR21

    ref_data = []
    with open(HERE / "data" / "gen_sigmas_sample_Si.xsections", "r") as f:
        for line in f:
            if not line.startswith("   ") or line.startswith("    "):
                continue
            ref_data.append([float(x) for x in line.strip().split()])

    ref_data = np.array(ref_data, dtype=KFLOAT)
    ref_data.shape = (20, 6)

    # Test regularized (small) log energy tables for Silicon,
    # For the values
    # Create medium here with same parameters, rho,sumA forced to match Mortran,
    # our elements data table has slightly different values (2.33, 28.0855)
    Si = Medium(1, "Si", rho=2.4, ap=0.01, up=50.0, mge=20, sumA=28.088)
    assert Si.gmfp01[0] == pytest.approx(ref_data[:, 0])
    assert Si.gmfp01[1] == pytest.approx(ref_data[:, 1])

    assert Si.gbr12[0, 0] == pytest.approx(ref_data[:, 2])
    assert Si.gbr12[0, 1] == pytest.approx(ref_data[:, 3])

    assert Si.gbr12[1, 0] == pytest.approx(ref_data[:, 4])
    assert Si.gbr12[1, 1] == pytest.approx(ref_data[:, 5])

