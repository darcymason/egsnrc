# media.py
"""Classes for Media and Elements"""
from dataclasses import dataclass
from enum import Enum
from math import exp, log
from typing import ClassVar
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
    ap: float
    "Lower photon cutoff energy"
    up: float
    rho: float = np.NaN
    "Medium density"
    cross_section_prefix: str = "xcom"
    "Upper photon cutoff energy"
    mge: int = MXGE
    "Number of energy intervals for interaction coefficients"
    photon_cross_sections: ClassVar[dict] = {}
    "Raw data from interaction cross-section files -dict[prefix][PHOTOELECTRIC/etc]"
    # Following are calculated on init
    # ge0: float = 0
    # ge1: float = 0
    # sigmas = {}  - dict with Interaction type key to array of sigmas

    def __post_init__(self):
        self.sigmas = {}
        if self.rho is np.NaN:
            if len(self.elements) == 1:
                self.rho = self.elements[0].rho
            else:
                raise NotImplementedError(
                    "Must supply a density rho for a multi-element Medium"
                )
        self.calc_all()

    def calc_all(self):
        self.calc_ge0_ge1()
        self.sumZ = sum(element.pz * element.z for element in self.elements)
        self.sumA = sum(element.pz * element.atomic_weight for element in self.elements)
        self.con1 = self.sumZ * self.rho / (self.sumA * 1.6605655)
        self.con2 = self.rho / (self.sumA * 1.6605655)
        self.calc_all_sigmas(self.cross_section_prefix)
        self.calc_branching_ratios()

    def calc_ge0_ge1(self):
        """Define log intervals over the lower and upper cutoff energies"""
        self.ge1 = (self.mge - 1) / log(self.up / self.ap)
        self.ge0 = 1 - self.ge1 * log(self.ap)

    def calc_all_sigmas(self, table_prefix):
        """Calculate all cross-sections for this medium

        Parameters
        ----------
        table_prefix : str
            Prefix of cross-section data file name, e.g. "xcom"
        """
        from egsnrc.hatch import DATA_DIR, get_xsection_table
        prefix = DATA_DIR / f"{table_prefix}_"
        intn = Interaction  # Just to shorten lines

        # Read data table files if not stored in the class already
        if table_prefix not in self.photon_cross_sections:
            tbls = self.photon_cross_sections[table_prefix] = {}
            tbls[intn.PHOTOELECTRIC] = get_xsection_table(f"{prefix}photo.data")
            tbls[intn.COMPTON] = get_xsection_table(f"{prefix}compton.data")
            tbls[intn.PAIR] = get_xsection_table(f"{prefix}pair.data")
            # XXX Rayleigh, Triplet, Photonuc

        # Get stored tables
        tbls = self.photon_cross_sections[table_prefix]
        photo_data = tbls[intn.PHOTOELECTRIC]
        compton_data = tbls[intn.COMPTON]
        pair_data = tbls[intn.PAIR]
        # XXX Rayleigh, Triplet, Photonuc

        self._calc_sigmas(intn.PHOTOELECTRIC, photo_data)
        self._calc_sigmas(intn.COMPTON, compton_data)
        self._calc_sigmas(intn.PAIR, pair_data)
        # XXX Rayleigh, Triplet, Photonuc

    @np.errstate(divide="raise")
    def _calc_sigmas(self, interaction, cross_sections):
        data = [0] * self.mge
        PAIR_OR_TRIPLET = (Interaction.PAIR, Interaction.TRIPLET)
        for element in self.elements:
            etmp, ftmp = cross_sections[element.z]

            # For pair or triplet, insert an extra data point for threshold energy
            # and change cross-sections (because ...?)
            if interaction in PAIR_OR_TRIPLET:
                eth = 2 * rm if interaction == Interaction.PAIR else 4 * rm
                ftmp = ftmp - 3 * np.log(1 - eth / np.exp(etmp))
                ftmp = np.insert(ftmp, 0, ftmp[0])
                etmp = np.insert(etmp, 0, log(eth))

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

        self.sigmas[interaction] = np.array(data)

    def calc_sigma(self, interaction, energy):
        gle = log(energy)
        sigmas = self.sigmas[interaction]
        fractional_index = (gle - log(self.ap)) * self.ge1
        index = int(fractional_index)
        p = fractional_index - index
        return (1 - p) * sigmas[index] + p * sigmas[index + 1]

    def calc_branching_ratios(self):
        list_2_to_mge = np.arange(2, self.mge + 1) # start 2 - used with diffs
        gle = (list_2_to_mge - self.ge0) / self.ge1
        #  exp_gle = np.exp(gle)
        # sig_KN = self.sumZ * egs_KN_sigma0(e)
        #     if ibcmp[1] > 1:
        #         if  input_compton_data:
                #     sig_KN = sig_compton(i)
                # else:
                #     #Apply the bound Compton correction to sig_KN"
                #     if  e <= bc_emin:
                #         bcf = exp(bc_data[1])
                #     elif  e < bc_emax:
                #         aj = 1 + log(e/bc_emin)/bc_dle
                #         j = int(aj)
                #         aj = aj - j
                #         bcf = exp(bc_data[j]*(1-aj) + bc_data[j+1]*aj)
                #     else:
                #         bcf = 1
                #     sig_KN = sig_KN*bcf
        # XXX for now, assume 'input_compton_data' for toy examples
        sig_KN = self.sigmas[Interaction.COMPTON]
        sig_pair = self.sigmas[Interaction.PAIR]
        sig_photo = self.sigmas[Interaction.PHOTOELECTRIC]

        sig_p  = sig_pair # XXX + sig_triplet
        sigma  = sig_KN + sig_p + sig_photo
        gmfp   = 1 / (sigma * self.con2)
        gbr1   = sig_p / sigma
        gbr2   = gbr1 + sig_KN / sigma
        # cohe   = sigma/(sig_rayleigh(i) + sigma)
        # photonuc = sigma/(sig_photonuc(i) + sigma)

        def arr1_0(arr):
            tmp1 = np.diff(arr) * self.ge1
            tmp0 = arr[1:] - tmp1 * gle
            result1 = np.append(tmp1, tmp1[-1]) # gmfp1(nge,med)=gmfp1(nge-1,med)
            result0 = np.append(tmp0, tmp0[-1])
            return result1, result0

        self.gmfp1, self.gmfp0 = arr1_0(gmfp)
        self.gbr11, self.gbr10 = arr1_0(gbr1)
        self.gbr21, self.gbr20 = arr1_0(gbr2)
        # XXX repeat for cohe and photonuc


if __name__ == "__main__":
    from egsnrc.hatch import DATA_DIR, get_xsection_table

    # PHOTO = Interaction.PHOTOELECTRIC
    # photo_data = get_xsection_table(DATA_DIR / "xcom_photo.data")
    # element_ta = Element(z=73, pz=1, wa=0)
    element_C = Element(z=6, pz=1)
    # c_en = photo_data[6][0]

    # # Try to match the energies of the original cross-section data
    # medium = Medium(
    #     "C", [element_C], rho=1, ap=exp(c_en[0]), up=exp(c_en[-1]), mge=len(c_en)
    # )
    medium = Medium("C", [element_C], ap=0.001, up=30.0)
    # sigmas = medium.sigmas[PHOTO]
    print(medium)

    # medium = Medium("Medium Ca", [Element(20, 1)], ap=0.001, up=2.0)
    for z in range(1, 100):
        Medium("test", [Element(z=z, pz=1)], ap=0.001, up=50.0)
    print(medium)

