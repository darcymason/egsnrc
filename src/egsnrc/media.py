# media.py
"""Classes for Media and Elements"""
from collections import namedtuple
from dataclasses import dataclass
from enum import Enum
from math import exp, log
from egsnrc.config import KINT, KFLOAT
from egsnrc.hatch import DATA_DIR, get_xsection_table
from egsnrc.elements import element_data, fcoulc, xsif
from egsnrc.params import MXGE
import numpy as np

from egsnrc.constants import REST_MASS


GBR_PAIR, GBR_COMPTON = np.arange(2, dtype=np.int32)

# Internal _Medium object for gpu calls, due to limitations of types that can be passed
_Medium = namedtuple("_Medium", "number rho gmfp01 ge gbr12")  # formula name


Vacuum = _Medium(
    np.int32(0),
    np.float32(0.0),  # "Vacuum", "Vacuum",
    np.empty((0, 0), dtype=np.float32),
    np.empty(0, dtype=np.float32),
    np.empty((0, 0, 0), dtype=np.float32),
)


Interaction = Enum(
    "Interaction",
    ["PHOTOELECTRIC", "RAYLEIGH", "COMPTON", "PAIR", "TRIPLET", "PHOTONUCLEAR"],
)


def egs_KN_sigma0(e_arr):
    """Calculate portion of Klein-Nishina sigma

    e_arr:  numpy.array
        Array of energies to convert
    """
    con = 0.1274783851
    ko_arr = e_arr / REST_MASS
    # XXX can rewrite to be more efficient by using numpy bool mask rather than looping
    result = []
    for ko, e in zip(ko_arr, e_arr):
        if ko < 0.01:
            # egs_KN_sigma0 = 8.*con/3.*(1-ko*(2-ko*(5.2-13.3*ko)))/prm
            result.append(
                8.0 * con / 3.0 * (1 - ko * (2 - ko * (5.2 - 13.3 * ko))) / REST_MASS
            )
        else:
            # c1 = 1./(ko*ko); c2 = 1. - 2*(1+ko)*c1; c3 = (1+2*ko)*c1
            # eps2 = 1; eps1 = 1./(1+2*ko)
            # egs_KN_sigma0 = (c1*(1./eps1-1./eps2) + c2*log(eps2/eps1)
            # + eps2*(c3+0.5*eps2)- eps1*(c3+0.5*eps1))/e*con
            c1 = 1.0 / (ko * ko)
            c2 = 1.0 - 2 * (1 + ko) * c1
            c3 = (1 + 2 * ko) * c1
            eps2 = 1.0
            eps1 = 1.0 / (1.0 + 2.0 * ko)
            result.append(
                (
                    c1 * (1.0 / eps1 - 1.0 / eps2)
                    + c2 * log(eps2 / eps1)
                    + eps2 * (c3 + 0.5 * eps2)
                    - eps1 * (c3 + 0.5 * eps1)
                )
                / e
                * con
            )
    return np.array(result)


@dataclass
class Medium:
    """Contain information for a specific physical medium"""

    number: int
    formula: str
    name: str = None
    rho: float = None
    "Medium density"
    ap: float = 0.001
    "Lower photon cutoff energy"
    up: float = 50.0
    "Upper photon cutoff energy"
    mge: int = MXGE
    "Number of energy intervals for interaction coefficients"
    table_prefix: str = "xcom"
    "Prefix for name of cross-section data files"
    sumA: float = None
    "Sum of atomic weights; mainly for test suite to match different element data"
    # Following are calculated on init
    # ge0: float = 0
    # ge1: float = 0
    # sigmas = {}  - dict with Interaction type key to array of sigmas

    def __post_init__(self):
        if self.name is None:
            self.name = self.formula
        if self.formula not in element_data:
            raise NotImplementedError("Currently only accept single chemical elements")
            # XXX can add library like chemparse, chemformula, or chempy to give us
            #    the proportions
        symbols = [self.formula]
        self.elements = [element_data[symbol] for symbol in symbols]
        self.proportions = [1]

        if self.rho is None:
            if len(self.elements) == 1:
                self.rho = self.elements[0].density
            else:
                raise NotImplementedError(
                    "Must supply a density rho for a multi-element Medium"
                )

        # Define log intervals over the lower and upper cutoff energies"""
        self.ge1 = (self.mge - 1) / log(self.up / self.ap)
        self.ge0 = 1 - self.ge1 * log(self.ap)
        self.ge = np.array(
            (self.ge0, self.ge1), KFLOAT
        )  # "Kernelize" here, user can inspect

        self.sumZ = sum(
            proportion * element.atomic_number
            for element, proportion in zip(self.elements, self.proportions)
        )
        if self.sumA is None:
            self.sumA = sum(
                proportions * element.atomic_weight
                for element, proportions in zip(self.elements, self.proportions)
            )
        # con1 = sumZ * rho / (self.sumA * 1.6605655)  # con1 never used?
        con2 = self.rho / (self.sumA * 1.6605655)

        # Calculate all cross-sections for this medium on regularized log lookup
        # First get physical data tables
        prefix = DATA_DIR / f"{self.table_prefix}_"
        intn = Interaction  # Just to shorten lines
        photo_data = get_xsection_table(f"{prefix}photo.data")
        compton_data = get_xsection_table(f"{prefix}compton.data")
        pair_data = get_xsection_table(f"{prefix}pair.data")
        triplet_data = get_xsection_table(f"{prefix}triplet.data")
        # XXX Rayleigh, Triplet, Photonuc

        self.sigmas = {}
        for interaction, data in (
            (intn.PHOTOELECTRIC, photo_data),
            (intn.COMPTON, compton_data),
            (intn.PAIR, pair_data),
            (intn.TRIPLET, triplet_data),
        ):  # XXX add Rayleigh, Triplet, Photonuc
            self.sigmas[interaction] = _calc_sigmas(
                interaction,
                data,
                self.elements,
                self.proportions,
                self.ge0,
                self.ge1,
                self.mge,
            )

        # Calculate branching ratios
        list_1_to_mge = np.arange(1, self.mge + 1)  # start 1 - used with diffs
        gle = (list_1_to_mge - self.ge0) / self.ge1
        self.energies = np.exp(gle)
        sig_KN = self.sumZ * egs_KN_sigma0(self.energies)
        #     if ibcmp[1] > 1:
        #         if  input_compton_data:
        #           sig_KN = sig_compton(i)
        #         else:
        #             # Apply the bound Compton correction to sig_KN"
        #             if  e <= bc_emin:
        #                bcf = exp(bc_data[1])
        #             elif  e < bc_emax:
        #                aj = 1 + log(e/bc_emin)/bc_dle
        #                j = int(aj)
        #                aj = aj - j
        #                bcf = exp(bc_data[j]*(1-aj) + bc_data[j+1]*aj)
        #             else:
        #                bcf = 1
        #             sig_KN = sig_KN*bcf
        # XXX for now, assume 'input_compton_data'
        # sig_KN = self.sigmas[Interaction.COMPTON]
        sig_pair = self.sigmas[Interaction.PAIR]
        sig_photo = self.sigmas[Interaction.PHOTOELECTRIC]
        sig_triplet = self.sigmas[Interaction.TRIPLET]

        sig_p = sig_pair + sig_triplet
        sigma = sig_KN + sig_p + sig_photo
        gmfp = 1 / (sigma * con2)
        gbr1 = sig_p / sigma
        gbr2 = gbr1 + sig_KN / sigma
        # cohe   = sigma/(sig_rayleigh(i) + sigma)
        # photonuc = sigma/(sig_photonuc(i) + sigma)

        def arr0_1(arr):
            tmp1 = np.diff(arr) * self.ge1
            tmp0 = arr[1:] - tmp1 * gle[1:]
            result1 = np.append(tmp1, tmp1[-1])  # gmfp1(nge,med)=gmfp1(nge-1,med)
            result0 = np.append(tmp0, tmp0[-1])
            return np.stack((result0, result1))

        # "Kernelize" these arrays here, so user can inspect actual values used
        self.gmfp01 = arr0_1(gmfp).astype(KFLOAT)
        gbr1 = arr0_1(gbr1).astype(KFLOAT)
        gbr2 = arr0_1(gbr2).astype(KFLOAT)
        self.gbr12 = np.stack((gbr1, gbr2))  # pair, compton, ...
        # XXX repeat for cohe and photonuc

    def kernelize(self):
        return _Medium(
            KINT(self.number),
            # name, formula,
            KFLOAT(self.rho),
            self.gmfp01,
            self.ge,
            self.gbr12,
        )

    # Not used except in test suite
    def calc_sigma(self, interaction, energy):
        gle = log(energy)
        sigmas = self.sigmas[interaction]
        fractional_index = (gle - log(self.ap)) * self.ge1
        index = int(fractional_index)
        p = fractional_index - index
        return (1 - p) * sigmas[index] + p * sigmas[index + 1]

    def calc_brems_pair(self):
        """Calculate factors related to brems and pair"""
        # Adapted from subroutine fix_brems in egsnrc.mortran
        Zt = 0
        Zb = 0
        Zf = 0
        for elem, pi in zip(self.elements, self.proportions):
            Zi = elem.atomic_number
            aux = pi * Zi * (Zi + xsif(Zi))
            Zt += aux
            Zb -= aux * log(Zi) / 3
            Zf += aux * fcoulc(Zi)

        Zv = (Zb - Zf) / Zt
        Zg = Zb / Zt
        fmax1 = 2 * (20.863 + 4 * Zg) - 2 * (20.029 + 4 * Zg) / 3
        fmax2 = 2 * (20.863 + 4 * Zv) - 2 * (20.029 + 4 * Zv) / 3

        # Initialize empty arrays for DL factors
        # Unlike in mortran code, make a 2d array to keep Medium tuple manageable
        # indexes 1 - 8, here keep an unused 0 in the 0th index
        dl = np.zeros(7, 9)
        bpar = np.zeros(3)

        # These are used in BREMS
        dl[1, 1] = (20.863 + 4 * Zg) / fmax1
        dl[2, 1] = -3.242 / fmax1
        dl[3, 1] = 0.625 / fmax1
        dl[4, 1] = (21.12 + 4 * Zg) / fmax1
        dl[5, 1] = -4.184 / fmax1
        dl[6, 1] = 0.952
        dl[1, 2] = (20.029 + 4 * Zg) / fmax1
        dl[2, 2] = -1.93 / fmax1
        dl[3, 2] = -0.086 / fmax1
        dl[4, 2] = (21.12 + 4 * Zg) / fmax1
        dl[5, 2] = -4.184 / fmax1
        dl[6, 2] = 0.952
        dl[1, 3] = (20.863 + 4 * Zv) / fmax2
        dl[2, 3] = -3.242 / fmax2
        dl[3, 3] = 0.625 / fmax2
        dl[4, 3] = (21.12 + 4 * Zv) / fmax2
        dl[5, 3] = -4.184 / fmax2
        dl[6, 3] = 0.952
        dl[1, 4] = (20.029 + 4 * Zv) / fmax2
        dl[2, 4] = -1.93 / fmax2
        dl[3, 4] = -0.086 / fmax2
        dl[4, 4] = (21.12 + 4 * Zv) / fmax2
        dl[5, 4] = -4.184 / fmax2
        dl[6, 4] = 0.952

        # and these in PAIR
        dl[1, 5] = 3 * (20.863 + 4 * Zg) - (20.029 + 4 * Zg)
        dl[2, 5] = 3 * (-3.242) - (-1.930)
        dl[3, 5] = 3 * (0.625) - (-0.086)
        dl[4, 5] = 2 * 21.12 + 8 * Zg
        dl[5, 5] = 2 * (-4.184)
        dl[6, 5] = 0.952
        dl[1, 6] = 3 * (20.863 + 4 * Zg) + (20.029 + 4 * Zg)
        dl[2, 6] = 3 * (-3.242) + (-1.930)
        dl[3, 6] = 3 * 0.625 + (-0.086)
        dl[4, 6] = 4 * 21.12 + 16 * Zg
        dl[5, 6] = 4 * (-4.184)
        dl[6, 6] = 0.952
        dl[1, 7] = 3 * (20.863 + 4 * Zv) - (20.029 + 4 * Zv)
        dl[2, 7] = 3 * (-3.242) - (-1.930)
        dl[3, 7] = 3 * (0.625) - (-0.086)
        dl[4, 7] = 2 * 21.12 + 8 * Zv
        dl[5, 7] = 2 * (-4.184)
        dl[6, 7] = 0.952
        dl[1, 8] = 3 * (20.863 + 4 * Zv) + (20.029 + 4 * Zv)
        dl[2, 8] = 3 * (-3.242) + (-1.930)
        dl[3, 8] = 3 * 0.625 + (-0.086)
        dl[4, 8] = 4 * 21.12 + 16 * Zv
        dl[5, 8] = 4 * (-4.184)
        dl[6, 8] = 0.952

        bpar[2] = dl[1, 7] / (3 * dl[1, 8] + dl[1, 7])
        bpar[1] = 12 * dl[1, 8] / (3 * dl[1, 8] + dl[1, 7])

        # Store arrays in the medium
        self.dl = dl
        self.bpar = bpar

        # DELCM = 136*m*exp(Zg), eq. (2.7.51)


@np.errstate(divide="raise")
def _calc_sigmas(interaction, cross_sections, elements, proportions, ge0, ge1, mge):
    data = [0] * mge
    PAIR_OR_TRIPLET = (Interaction.PAIR, Interaction.TRIPLET)
    for element, proportion in zip(elements, proportions):
        etmp, ftmp = cross_sections[element.atomic_number]

        # For pair or triplet, insert an extra data point for threshold energy
        # and change cross-sections (because ...?)
        if interaction in PAIR_OR_TRIPLET:
            eth = 2 * REST_MASS if interaction == Interaction.PAIR else 4 * REST_MASS
            ftmp = ftmp - 3 * np.log(1 - eth / np.exp(etmp))
            ftmp = np.insert(ftmp, 0, ftmp[0])
            etmp = np.insert(etmp, 0, log(eth))

        for k in range(mge):
            gle = (k + 1 - ge0) / ge1
            energy = exp(gle)

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
                    p = (energy - exp(etmp[kk])) / (exp(etmp[kk + 1]) - exp(etmp[kk]))
                    sig = p * exp(ftmp[kk + 1]) + (1 - p) * exp(ftmp[kk])

            if interaction in PAIR_OR_TRIPLET and energy > eth:
                sig = sig * (1 - eth / energy) ** 3

            data[k] += proportion * sig

    return np.array(data, dtype=KFLOAT)
