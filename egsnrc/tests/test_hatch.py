import pytest
from pathlib import Path
from egsnrc.hatch import read_pegs


HERE  = Path(__file__).resolve().parent
TEST_DATA = HERE / "data"


def test_pegs_read():
    """Correct values are returned on reading PEGS4 data file"""
    requested_media = {'TA', 'NaI', 'H2O'}
    media = read_pegs(TEST_DATA / "tutor_data.pegs4dat", requested_media)

    nai = media['nai']
    assert nai['medium'] == "NAI"
    assert media.keys() == set(x.lower() for x in requested_media)
    assert 199 == nai['mge']
    assert 46 == nai['meke']

    # brempr
    # DL1-6, first  0.999978E+00 -0.251878E+00  0.578002E-01  0.995067E+00 -0.270359E+00  0.952000E+00
    # DL1-6, last  0.100044E+01 -0.210042E+00  0.322380E-01  0.103054E+01 -0.301584E+00  0.952000E+00
    # DL1-6, second last  0.100028E+01 -0.230911E+00  0.445155E-01  0.101830E+01 -0.298004E+00  0.952000E+00

    # Values exported from tutor4 Mortran via $egs_info macro for TA:
    ta = media['ta']
    tol = 1e-5
    # Sample expected values:
    assert 0.999978E+00 == pytest.approx(ta['dl1'][0], tol)
    assert 0.100044E+01 == pytest.approx(ta['dl1'][-1], tol)
    assert 0.952000E+00 == pytest.approx(ta['dl6'][-1], tol)
    assert 0.445155E-01 == pytest.approx(ta['dl3'][-2], tol)
    assert 0.952000E+00 == pytest.approx(ta['dl6'][-2], tol)

    # elecin
    # esig0(1,1), psig0(1,1),ededx0(1,1), pdedx0(1,1), EBR10, PBR10, PBR20, TMXS0
    # and similar for esig1 etc, and then all again with (neke,1) indexes

    # 0.671753E+01    0.954806E+01    0.307003E+01    0.187255E+01    0.100000E+01    0.566277E+00    0.104216E+01    0.180656E-02
    # 0.524882E+01    0.500037E+01    0.640958E+01    0.639328E+01    0.722474E+00    0.740728E+00    0.992775E+00   -0.226768E+02
    # 0.123095E+01    0.258473E+01   -0.425757E+01   -0.512867E+01    0.000000E+00   -0.194480E+00    0.911563E-01    0.800434E-03
    # 0.192685E+01    0.197515E+01    0.182861E+00    0.186221E+00    0.263692E-01    0.226901E-01    0.147533E-02    0.688693E+01

    var_names = "esig psig ededx pdedx ebr1 pbr1 pbr2 tmxs".split()
    vars_0 = [var_name + "0" for var_name in var_names]
    vars_1 = [var_name + "1" for var_name in var_names]
    expected_var0_0 = (
        "0.671753E+01  0.954806E+01  0.307003E+01  0.187255E+01 "
        "0.100000E+01  0.566277E+00  0.104216E+01  0.180656E-02"
    ).split()
    expected_var0_last = (
        "0.524882E+01    0.500037E+01    0.640958E+01    0.639328E+01 "
        "0.722474E+00    0.740728E+00    0.992775E+00   -0.226768E+02"
    ).split()

    expected_var1_0 = (
        "0.123095E+01    0.258473E+01   -0.425757E+01   -0.512867E+01 "
        "0.000000E+00   -0.194480E+00    0.911563E-01    0.800434E-03"
    ).split()

    expected_var1_last = (
        "0.192685E+01    0.197515E+01    0.182861E+00    0.186221E+00 "
        "0.263692E-01    0.226901E-01    0.147533E-02    0.688693E+01"
    ).split()

    for var, expected in zip(vars_0, expected_var0_0):
        assert float(expected) == pytest.approx(ta[var][0])

    for var, expected in zip(vars_0, expected_var0_last):
        assert float(expected) == pytest.approx(ta[var][-1])

    for var, expected in zip(vars_1, expected_var1_0):
        assert float(expected) == pytest.approx(ta[var][0])

    for var, expected in zip(vars_1, expected_var1_last):
        assert float(expected) == pytest.approx(ta[var][-1])

    # PHOTIN
    # EBINDA, GE0, GE1    0.674160E-01    0.108496E+03    0.231254E+02
    assert 0.674160E-01 == ta['ebinda']
    assert 0.108496E+03 == ta['ge0']
    assert 0.231254E+02 == ta['ge1']

    var_names = 'gmfp gbr1 gbr2'.split()
    vars_0 = [var_name + "0" for var_name in var_names]
    vars_1 = [var_name + "1" for var_name in var_names]
    expected_var0_0 =  "0.924862E-02    0.000000E+00    0.968388E-02".split()
    expected_var0_last = "0.276774E+01    0.774469E+00    0.995661E+00".split()
    expected_var1_0 =  "0.187401E-02    0.000000E+00    0.196035E-02".split()
    expected_var1_last = "-0.205889E+00    0.453419E-01    0.904219E-03".split()
    for var, expected in zip(vars_0, expected_var0_0):
        assert float(expected) == pytest.approx(ta[var][0])

    for var, expected in zip(vars_0, expected_var0_last):
        assert float(expected) == pytest.approx(ta[var][-1])

    for var, expected in zip(vars_1, expected_var1_0):
        assert float(expected) == pytest.approx(ta[var][0])

    for var, expected in zip(vars_1, expected_var1_last):
        assert float(expected) == pytest.approx(ta[var][-1])

