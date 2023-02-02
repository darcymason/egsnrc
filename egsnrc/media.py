# media.py
"""Classes for Media and Elements"""
from collections import namedtuple
from dataclasses import dataclass
from enum import Enum
from math import exp, log
from typing import ClassVar
from egsnrc.hatch import DATA_DIR, get_xsection_table
from egsnrc.elements import element_data
from egsnrc.params import MXGE
import numba as nb

import numpy as np

from egsnrc.constants import REST_MASS


GBR_PAIR, GBR_COMPTON = np.arange(2, dtype=np.int32)

# Internal _Medium object for gpu calls, due to limitations of types that can be passed
Medium = namedtuple("_Medium", "number formula name rho gmfp ge gbr")


Vacuum = Medium(
    np.int32(0), "Vacuum", "", np.float32(0.0),
    np.array([], dtype=np.float32),
    np.array([], dtype=np.float32),
    np.array([], dtype=np.float32),
)


def _media_to_arrays(medium_list):
    """Given a Python list of Medium class instances, return a tuple of _medium_obj's"""
    #  This is used internally to reformat calls into GPU code, due to current
    #     limitations on types that can be passed
    # Media numbers should normally be close together, but don't assume that
    # Allocate enough array for the highest number.
    medium_numbers = [medium.medium_number for medium in medium_list]
    medium_objs = [_medium_obj(medium) for medium in medium_list]

    # Create 'sparse' tuple so can index directly by medium number
    # XXX could be better perhaps to use item 0 to store indices to redirect
    media = tuple(
        medium_objs[i] if i in medium_numbers else 0
        for i in range(max(medium_numbers) + 1)
    )
    return media

def _medium_obj(medium):
    """Return a namedtuple containing the Python Medium class information"""
    m = medium
    _float = nb.float32
    # _int = nb.int32
    return _Medium(_float(m.rho), m.ge, m.gbr)


Interaction = Enum(
    "Interaction",
    ["PHOTOELECTRIC", "RAYLEIGH", "COMPTON", "PAIR", "TRIPLET", "PHOTONUCLEAR"],
)


def make_medium(
    number, formula, name=None, rho=None, ap=0.001, up=50.0, mge=MXGE,
    table_prefix="xcom"
):
    """Given a chemical formula*, return a Medium object with the interaction data

    *Currently only takes a single chemical element
    """
    if name is None:
        name = formula
    if not formula in element_data:
        raise NotImplementedError("Currently only accept single chemical elements")
        # XXX can add library like chemparse, chemformula, or chempy to give us
        #    the proportions
    symbols = [formula]
    elements = [element_data[symbol] for symbol in symbols]
    proportions = [1]
    sigmas = {}
    if rho is None:
        if len(elements) == 1:
            rho = elements[0].density
        else:
            raise NotImplementedError(
                "Must supply a density rho for a multi-element Medium"
            )

    # Define log intervals over the lower and upper cutoff energies"""
    ge1 = (mge - 1) / log(up / ap)
    ge0 = 1 - ge1 * log(ap)
    ge = np.array((ge0, ge1), dtype=np.float32)

    # The following (except con2) aren't used at the moment?
    sumZ = sum(
        proportion * element.atomic_number
        for element, proportion in zip(elements, proportions)
    )
    sumA = sum(
        proportions * element.atomic_weight
        for element, proportions in zip(elements, proportions)
    )
    con1 = sumZ * rho / (sumA * 1.6605655)  # con1 never used?
    con2 = rho / (sumA * 1.6605655)

    # Calculate all cross-sections for this medium on regularized log lookup
    # First get physical data tables
    prefix = DATA_DIR / f"{table_prefix}_"
    intn = Interaction  # Just to shorten lines
    photo_data = get_xsection_table(f"{prefix}photo.data")
    compton_data = get_xsection_table(f"{prefix}compton.data")
    pair_data = get_xsection_table(f"{prefix}pair.data")
    # XXX Rayleigh, Triplet, Photonuc

    sigmas = {}
    for interaction, data in (
        (intn.PHOTOELECTRIC, photo_data), (intn.COMPTON, compton_data),
        (intn.PAIR, pair_data)
    ):  # XXX add Rayleigh, Triplet, Photonuc
        sigmas[interaction] = _calc_sigmas(
            interaction, data, elements, proportions, ge0, ge1, mge
        )

    # Calculate branching ratios
    list_2_to_mge = np.arange(2, mge + 1) # start 2 - used with diffs
    gle = (list_2_to_mge - ge0) / ge1
    #  exp_gle = np.exp(gle)
    # sig_KN = sumZ * egs_KN_sigma0(e)
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
    # XXX for now, assume 'input_compton_data'
    sig_KN = sigmas[Interaction.COMPTON]
    sig_pair = sigmas[Interaction.PAIR]
    sig_photo = sigmas[Interaction.PHOTOELECTRIC]

    sig_p  = sig_pair # XXX + sig_triplet
    sigma  = sig_KN + sig_p + sig_photo
    gmfp   = 1 / (sigma * con2)
    gbr1   = sig_p / sigma
    gbr2   = gbr1 + sig_KN / sigma
    # cohe   = sigma/(sig_rayleigh(i) + sigma)
    # photonuc = sigma/(sig_photonuc(i) + sigma)

    def arr0_1(arr):
        tmp1 = np.diff(arr) * ge1
        tmp0 = arr[1:] - tmp1 * gle
        result1 = np.append(tmp1, tmp1[-1]) # gmfp1(nge,med)=gmfp1(nge-1,med)
        result0 = np.append(tmp0, tmp0[-1])
        return np.stack((result0, result1))

    gmfp = arr0_1(gmfp)
    gbr1 = arr0_1(gbr1)
    gbr2 = arr0_1(gbr2)
    gbr = np.stack((gbr1, gbr2))  # pair, compton, ...
    # XXX repeat for cohe and photonuc

    return Medium(number, name, formula, rho, gmfp, ge, gbr)


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

            data[k] += proportion * sig

    return np.array(data)

# def calc_sigma(sigmas, interaction, energy, ap, ge1):
#     gle = log(energy)
#     sigmas = sigmas[interaction]
#     fractional_index = (gle - log(ap)) * ge1
#     index = int(fractional_index)
#     p = fractional_index - index
#     return (1 - p) * sigmas[index] + p * sigmas[index + 1]



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

