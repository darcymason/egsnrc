# elements.py
"""Information about Periodic Table elements
"""
#  From pegs_page.cpp 2023-01-03

from math import log
from typing import NamedTuple


class ElemData(NamedTuple):
    atomic_number: int
    atomic_weight: float
    mean_excitation_energy: float
    density: float

class Zdata(NamedTuple):
    symbol: str
    atomic_weight: float
    mean_excitation_energy: float
    density: float


element_data = {
    "H": ElemData(1, 1.0079, 19.2, 8.3748e-05),
    "He": ElemData(2, 4.0026, 41.8, 0.00016632),
    "Li": ElemData(3, 6.941, 40, 0.534),
    "Be": ElemData(4, 9.0122, 63.7, 1.848),
    "B": ElemData(5, 10.812, 76, 2.37),
    "C": ElemData(6, 12.011, 78, 2),
    "N": ElemData(7, 14.0067, 82, 0.0011653),
    "O": ElemData(8, 15.9994, 95, 0.0013315),
    "F": ElemData(9, 18.9984, 115, 0.0015803),
    "Ne": ElemData(10, 20.1797, 137, 0.0008385),
    "Na": ElemData(11, 22.9898, 149, 0.971),
    "Mg": ElemData(12, 24.305, 156, 1.74),
    "Al": ElemData(13, 26.9815, 166, 2.702),
    "Si": ElemData(14, 28.0855, 173, 2.33),
    "P": ElemData(15, 30.9738, 173, 2.2),
    "S": ElemData(16, 32.066, 180, 2),
    "Cl": ElemData(17, 35.4527, 174, 0.0029947),
    "Ar": ElemData(18, 39.948, 188, 0.001662),
    "K": ElemData(19, 39.0983, 190, 0.862),
    "Ca": ElemData(20, 40.078, 191, 1.55),
    "Sc": ElemData(21, 44.9559, 216, 2.989),
    "Ti": ElemData(22, 47.88, 233, 4.54),
    "V": ElemData(23, 50.9415, 245, 6.11),
    "Cr": ElemData(24, 51.9961, 257, 7.18),
    "Mn": ElemData(25, 54.938, 272, 7.44),
    "Fe": ElemData(26, 55.847, 286, 7.874),
    "Co": ElemData(27, 58.9332, 297, 8.9),
    "Ni": ElemData(28, 58.69, 311, 8.902),
    "Cu": ElemData(29, 63.546, 322, 8.96),
    "Zn": ElemData(30, 65.39, 330, 7.133),
    "Ga": ElemData(31, 69.723, 334, 5.904),
    "Ge": ElemData(32, 72.61, 350, 5.323),
    "As": ElemData(33, 74.9216, 347, 5.73),
    "Se": ElemData(34, 78.96, 348, 4.5),
    "Br": ElemData(35, 79.904, 357, 0.0070722),
    "Kr": ElemData(36, 83.8, 352, 0.0034783),
    "Rb": ElemData(37, 85.4678, 363, 1.532),
    "Sr": ElemData(38, 87.62, 366, 2.54),
    "Y": ElemData(39, 88.9059, 379, 4.469),
    "Zr": ElemData(40, 91.224, 393, 6.506),
    "Nb": ElemData(41, 92.9064, 417, 8.57),
    "Mo": ElemData(42, 95.94, 424, 10.22),
    "Tc": ElemData(43, 97.9072, 428, 11.5),
    "Ru": ElemData(44, 101.07, 441, 12.41),
    "Rh": ElemData(45, 102.906, 449, 12.41),
    "Pd": ElemData(46, 106.42, 470, 12.02),
    "Ag": ElemData(47, 107.868, 470, 10.5),
    "Cd": ElemData(48, 112.411, 469, 8.65),
    "In": ElemData(49, 114.82, 488, 7.31),
    "Sn": ElemData(50, 118.71, 488, 7.31),
    "Sb": ElemData(51, 121.75, 487, 6.691),
    "Te": ElemData(52, 127.6, 485, 6.24),
    "I": ElemData(53, 126.904, 491, 4.93),
    "Xe": ElemData(54, 131.29, 482, 0.0054854),
    "Cs": ElemData(55, 132.905, 488, 1.873),
    "Ba": ElemData(56, 137.327, 491, 3.5),
    "La": ElemData(57, 138.905, 501, 6.154),
    "Ce": ElemData(58, 140.115, 523, 6.657),
    "Pr": ElemData(59, 140.908, 535, 6.71),
    "Nd": ElemData(60, 144.24, 546, 6.9),
    "Pm": ElemData(61, 144.913, 560, 7.22),
    "Sm": ElemData(62, 150.36, 574, 7.46),
    "Eu": ElemData(63, 151.965, 580, 5.243),
    "Gd": ElemData(64, 157.25, 591, 7.9004),
    "Tb": ElemData(65, 158.925, 614, 8.229),
    "Dy": ElemData(66, 162.5, 628, 8.55),
    "Ho": ElemData(67, 164.93, 650, 8.795),
    "Er": ElemData(68, 167.26, 658, 9.066),
    "Tm": ElemData(69, 168.934, 674, 9.321),
    "Yb": ElemData(70, 173.04, 684, 6.73),
    "Lu": ElemData(71, 174.967, 694, 9.84),
    "Hf": ElemData(72, 178.49, 705, 13.31),
    "Ta": ElemData(73, 180.948, 718, 16.654),
    "W": ElemData(74, 183.85, 727, 19.3),
    "Re": ElemData(75, 186.207, 736, 21.02),
    "Os": ElemData(76, 190.2, 746, 22.57),
    "Ir": ElemData(77, 192.22, 757, 22.42),
    "Pt": ElemData(78, 195.08, 790, 21.45),
    "Au": ElemData(79, 196.966, 790, 19.32),
    "Hg": ElemData(80, 200.59, 800, 13.546),
    "Tl": ElemData(81, 204.383, 810, 11.72),
    "Pb": ElemData(82, 207.2, 823, 11.34),
    "Bi": ElemData(83, 208.98, 823, 9.747),
    "Po": ElemData(84, 208.982, 830, 9.32),
    "At": ElemData(85, 209.987, 825, 9.32),
    "Rn": ElemData(86, 222.018, 794, 0.0090662),
    "Fr": ElemData(87, 223.02, 827, 1),
    "Ra": ElemData(88, 226.025, 826, 5),
    "Ac": ElemData(89, 227.028, 841, 10.07),
    "Th": ElemData(90, 232.038, 847, 11.72),
    "Pa": ElemData(91, 231.036, 878, 15.37),
    "U": ElemData(92, 238.029, 890, 18.95),
    "Np": ElemData(93, 237.048, 902, 20.25),
    "Pu": ElemData(94, 239.052, 921, 19.84),
    "Am": ElemData(95, 243.061, 934, 13.67),
    "Cm": ElemData(96, 247.07, 939, 13.51),
    "Bk": ElemData(97, 247.07, 952, 14),
    "Cf": ElemData(98, 251.08, 966, 10),
    "Es": ElemData(99, 252.083, 980, 10),
    "Fm": ElemData(100, 257.095, 994, 10),
}


# Reorganize dict by atomic number in case someone wants to use that way
element_data_by_z = {
    v.atomic_number: Zdata(k, v.atomic_weight, v.mean_excitation_energy, v.density)
    for k, v in element_data.items()
}


FINE = 137.03604
# Data for ALRAD and ALRADP for element z <= 4
# taken from Table B.2 in Y.Tsai Rev.Mod.Phys. 46,815(1974)
ALRAD = (05.31, 4.79, 4.74, 4.71)
ALRADP = (6.144, 5.621, 5.805, 5.924)
A1440 = 1194.0
A183 = 184.15


def fcoulc(Z):
    """Coulomb correction term in pair production and bremsstrahlung crosssections."""
    asq = Z / FINE
    asq = asq * asq
    return asq * (
        1.0 / (1.0 + asq) + 0.20206 + asq * (
            -0.0369 + asq * (0.0083 + asq * (-0.002))
        )
    )

def xsif(Z):
    """Function to account for brem and pair prodn in the field of atomic electrons."""
    if Z <= 4:
        iZ = int(Z) - 1  # -1 in Python for indexing
        return ALRADP(iZ) / (ALRAD(iZ) - fcoulc(Z))
    else:
        return log(A1440 * Z**(-0.666667)) / (log(A183 * Z**(-0.33333)) - fcoulc(Z))
