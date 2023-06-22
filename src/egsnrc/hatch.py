from pathlib import Path
import numpy as np

# hatch.py
"""Get physics data for EGSnrc run

Hatch is called before simulation begins.
For Photons, hatch calls egs_init_user_photon, which in turn
opens files via egsi_get_data, for compton, photoelectric, pair, triplet,
Rayleigh (depending on the settings) and does corrections on some of the data.
"""

from pathlib import Path
import numpy as np
import logging

logger = logging.getLogger("egsnrc")

DATA_DIR = Path(__file__).resolve().parent / "data"


def get_xsection_table(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    i_line = 0
    data = {}
    for z in range(1, 101):
        count = int(lines[i_line].split("#")[0])
        i_line += 1
        # 2 values per item, 8 values stored per line, so 4 data point pairs
        # Calc number of lines needed
        data_lines = count // 4 + (1 if count % 4 else 0)
        z_data = np.loadtxt(
            x for i in range(data_lines) for x in lines[i_line + i].strip().split()
        )
        # Reformat so have (array of energies, array of sigmas)
        z_data = z_data.reshape((-1, 2)).transpose()
        data[z] = z_data
        # print(f"Count {count}, len(data): {len(z_data)}")
        i_line += data_lines

    return data
