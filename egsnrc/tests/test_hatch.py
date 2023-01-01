import pytest
from egsnrc.hatch import DATA_DIR, photo_binding_energies, get_xsection_table


def test_binding_energies():
    """Test photo binding energies calculation from xcom photo data above 1 keV"""
    # Known binding energies (eV) from https://xdb.lbl.gov/Section1/Table_1-1.pdf
    known_binding_energies = {
        1: "13.6",
        6: "284.2",
        10: "870.2 48.5 21.7 21.6",
        20: "4038.5 438.4 349.7 346.2 44.3 25.4 25.4",
        73: (
            "67416 11682 11136 9881 2708 2469 2194 1793 1735 563.4 463.4 400.9 "
            "237.9 226.4 23.5 21.6 69.7 42.2 32.7"
        ),
        82: (
            "88005 15861 15200 13035 3851 3554 3066 2586 2484 891.8 761.9 643.5 "
            "434.3 412.2 141.7 136.9 147 106.4 83.3 20.7 18.1"
        ),
        92: (
            "115606 21757 20948 17166 5548 5182 4303 3728 3552 1439 1271 1043 "
            "778.3 736.2 388.2 377.4 321 257 192 102.8 94.2 43.9 26.8 16.8"
        ),
    }
    photo_data = get_xsection_table(DATA_DIR / "xcom_photo.data")
    binding_energies = photo_binding_energies(photo_data)
    for Z, known_energies_str in known_binding_energies.items():
        known_energies = sorted(
            float(x) / 1_000_000 for x in known_energies_str.split() if float(x) >= 1000
        )
        assert binding_energies[Z] == pytest.approx(known_energies, abs=0.00001)
