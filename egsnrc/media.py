# media.py
"""Classes for Media and Elements"""
from dataclasses import dataclass
from enum import Enum
from math import exp, log
from egsnrc.elements import element_data

import numpy as np

# from egsnrc.commons import rm
rm = 0.51099896430969238  # XXX need to put in global constants


MXGE = 2000  # EgsNRC default value; Number of energy intervals for sigma calcs


@dataclass
class Element:
    z: int
    """Atomic number"""
    pz: float
    """Proportion of element of type i.
       If a compound, the number of atoms of the element in the molecule.
       If a mixture,such as concrete, pz could be the percent of the atoms
    """

    # Following are filled in by standard tables if not supplied
    symbol: str = ""
    atomic_weight: float = 0
    rho: float = 0  # density

    def __post_init__(self):
        elem = element_data[self.z]
        if not self.symbol:
            self.symbol = elem.symbol
        if self.atomic_weight == 0:
            self.atomic_weight = elem.atomic_weight
        if self.rho == 0:
            self.rho = elem.density

Interaction = Enum(
    "Interaction",
    ["PHOTOELECTRIC", "RAYLEIGH", "COMPTON", "PAIR", "TRIPLET", "PHOTONUCLEAR"],
)


@dataclass  # XXX make frozen=True later - not post_init does calculations
class Medium:
    """Contain information for a specific physical medium"""

    name: str
    elements: list[Element]
    rho: float
    "Medium density"
    ap: float
    "Lower photon cutoff energy"
    up: float
    "Upper photon cutoff energy"
    mge: int = MXGE
    "Number of energy intervals for interaction coefficients"
    # Following are calculated on init
    # ge0: float = 0
    # ge1: float = 0
    # sigmas = {}
    # sig_photo: float = 0
    # sig_rayleigh: float = 0
    # sig_pair: float = 0
    # sig_triplet: float = 0
    # sig_photonuc: float = 0

    def __post_init__(self):
        self.sigmas = {}
        self.calc_all()

    def calc_all(self):
        self.calc_ge0_ge1()
        self.sumZ = sum(element.pz * element.z for element in self.elements)
        self.sumA = sum(element.pz * element.atomic_weight for element in self.elements)
        self.con1 = self.sumZ * self.rho / (self.sumA * 1.6605655)
        self.con2 = self.rho / (self.sumA * 1.6605655)


    def calc_ge0_ge1(self):
        """Define log intervals over the lower and upper cutoff energies"""
        self.ge1 = (self.mge - 1) / log(self.up / self.ap)
        self.ge0 = 1 - self.ge1 * log(self.ap)
        print(f"{self.mge=}  {self.ge1=}  {self.ge0=}")

    def calc_sigmas(self, interaction, cross_sections):
        data = [0] * self.mge
        PAIR_OR_TRIPLET = (Interaction.PAIR, Interaction.TRIPLET)
        for element in self.elements:
            etmp, ftmp = cross_sections[element.z]

            # For pair or triplet, insert an extra data point for threshold energy
            # and change cross-sections (because ...?)
            if interaction in PAIR_OR_TRIPLET:
                eth = 2 * rm if interaction == Interaction.PAIR else 4 * rm
                ftmp = ftmp - 3 * np.log(1 - eth / np.exp(etmp))
                ftmp.insert(0, ftmp[0])
                etmp.insert(0, log(eth))

            for k in range(self.mge):
                gle = (k + 1 - self.ge0) / self.ge1
                exp_gle = exp(gle)
                # Check within bounds
                if not etmp[0] <= gle < etmp[-1]:  # not within
                    if interaction in PAIR_OR_TRIPLET:
                        sig = 0 if gle < etmp[0] else exp(ftmp[-1])
                    elif interaction == Interaction.PHOTONUCLEAR:
                        sig = 0
                    else:
                        raise LookupError(
                            f"Energy {exp(gle)} is outside the available "
                            f"data range of {exp(etmp[0])} "
                            f"to {exp(etmp[-1])}'"
                        )
                else:  # within bounds of cross-section data
                    # searchsorted right:  a[i-1] <= v < a[i]
                    kk = etmp.searchsorted(gle, side="right") - 1
                    if interaction != Interaction.PHOTONUCLEAR:
                        # log/log interpolation
                        p = (gle - etmp[kk]) / (etmp[kk + 1] - etmp[kk])
                        sig = exp(p * ftmp[kk + 1] + (1 - p) * ftmp[kk])
                    else:
                        #  lin/lin interpolation for photonuc"
                        p = (exp_gle - exp(etmp[kk])) / (
                            exp(etmp[kk + 1]) - exp(etmp[kk])
                        )
                        sig = p * exp(ftmp[kk + 1]) + (1 - p) * exp(ftmp[kk])

                if interaction in PAIR_OR_TRIPLET and exp_gle > eth:
                    sig = sig * (1 - eth / exp_gle) ** 3

                data[k] += element.pz * sig

        self.sigmas[interaction] = data

    def calc_sigma(self, interaction, energy):
        gle = log(energy)
        sigmas = self.sigmas[interaction]
        fractional_index = (gle - log(self.ap)) * self.ge1
        index = int(fractional_index)
        p = fractional_index - index
        return (1 - p) * sigmas[index] + p * sigmas[index + 1]


if __name__ == "__main__":
    from egsnrc.hatch import DATA_DIR, get_xsection_table

    PHOTO = Interaction.PHOTOELECTRIC
    photo_data = get_xsection_table(DATA_DIR / "xcom_photo.data")
    # element_ta = Element(z=73, pz=1, wa=0)
    element_C = Element(z=6, pz=1, wa=0)
    c_en = photo_data[6][0]

    # Try to match the energies of the original cross-section data
    medium = Medium(
        "C", [element_C], rho=1, ap=exp(c_en[0]), up=exp(c_en[-1]), mge=len(c_en)
    )
    # medium = Medium("C", [element_C], rho=1, ap=0.001, up=50.0, mge=100)
    medium.calc_sigmas(PHOTO, photo_data)
    logE, log_sig = photo_data[6]
    sig = np.exp(log_sig)
    photo_E = np.exp(logE)
    sigmas = medium.sigmas[PHOTO]
    print(medium)

    print("  Log e: ", logE[:3], "...", logE[-3:])
    print("      e: ", photo_E[:3], "...", photo_E[-3:])
    print(" sigmas: ", sigmas[:3], "...", sigmas[-3:])
    print("tbl sig:", sig[:3], "...", sig[-3:])
    print(f"{len(sig)=}   {len(sigmas)=}")

