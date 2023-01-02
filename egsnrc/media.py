# media.py
"""Classes for Media and Elements"""
from dataclasses import dataclass
from enum import Enum
from math import exp, log

import numpy as np

# from egsnrc.commons import rm
rm = 0.51099896430969238


MXGE = 2000  # EgsNRC default value; Number of energy intervals for sigma calcs

# From Subroutine MIX:
# -------------------
#  NE       NUMBER OF DIFFERENT TYPES OF ATOMS IN THE MATERIAL.
#  PZ(I)    PROPORTION OF ELEMENT OF TYPE I.  IF A COMPOUND,
#           THEN PZ(I) WILL BE THE NUMBER OF ATOMS OF TYPE I IN THE MOLECULE.
#           IF A MIXTURE,SUCH AS CONCRETE, PZ(I) COULD BE THE PER CENT OF
#           THE ATOMS WHICH ARE OF TYPE I.
#  Z(I)     PERIODIC NUMBER OF ATOMS OF TYPE I
#  WA(I)    ATOMIC WEIGHT FOR ATOMS OF TYPE I.
#  WM = SUM(PZ(I)*WA(I)) = MOLECULAR  WEIGHT IF A COUMPOUND
#           OR A 'MIXTURE WEIGHT' IF A MIXTURE.
#  RHO      DENSITY OF THE MATERIAL. (IN GRAMS/CM**3)
#  RHOZ(I)  PARTIAL DENSITY DUE TO ATOMS OF TYPE I. (GM/CM**3)
#           ELECTRON DENSITY VARIABLE
#  ZC = SUM(PZ(I)*Z(I)) = NUMBER OF ELECTRONS/MOLECULE
#           BREMSSTRAHLUNG AND PAIR PRODUCTION VARIABLES ARE WEIGHTE
#  BY PZ(I)*Z(I)**2 FOR THE NUCLEUS, AND BY PZ(I)*Z(I)*XSI(I) FOR
#           ATOMIC ELECTRONS.
#  TPZ = SUM(PZ(I))
#  XSI(I) = LOG(A1440/Z(I)**(2./3.))/(LOG(A183/Z(I)**(1./3.))  -
#                FCOUL(Z(I)) )
#  ZZX(I) =  PZ(I)*Z(I)*(Z(I)+XSI(I)) = BREMS AND PAAR WEIGHTS
#  EZ = ZC/TPZ  EFFECTIVE Z
#  ZT = SUM(ZZX(I))
#  ZA = LOG(A183)*ZT   BUTCHER AND MESSELS L.C.'A' (1960)P.18
#  ZB = SUM(ZZX(I)*LOG(Z(I)**(-1./3.)  B&M'S L.C.'B' IBID.
#  ZF = SUM(ZZX(I)*FCOUL(Z(I))),WHERE FCOUL IS THE COULOMB
#           CORRECTION FUNCTION.
#  RATIOS--
#  ZG = ZB/ZT ,EXP(ZG)=WEIGHTED GEOMETRIC MEAN OF Z**(-1/3)
#  ZP = ZB/ZA , B&M IBID.P18 L.C.'P'
#  ZV= (ZB-ZF)/ZT
#  ZU = (ZB-ZF)/ZA
# ...
#  ZZ(I) = PZ(I)*Z(I)*(Z(I)+$FUDGEMS)
#  ZS = SUM(ZZ(I))
#  ZE = SUM(ZZ(I)*LOG(Z(I)**(-2./3.)))
#  ZX = SUM(ZZ(I)*LOG(1.+3.34*(FSC*Z(I))**2))
#                ELECTON DENSITY(ELECTRONS/CM**3)
#  EDEN=AN*RHO/WM*ZC
#           RADIATION LENGTH
#  USEFUL FOR GAUGING THE STEP SIZE, EVEN IF IT IS NOT USED AS THE
#  UNIT OF DISTANCE.
#   1./RLC =(AN*RHO/WM)*4.0*FSC*R0**2*
#     SUM( Z(I)*(Z(I)+XSI(I))*(LOG(A183*Z(I)**(-1./3.)-FCOUL(Z(I)))
#         =(AN*RHO/WM)*4.*FSC*R0**2*(ZAB-ZF)


@dataclass
class Element:
    z: int
    """Atomic number"""
    pz: float
    """Proportion of element of type i.
       If a compound, the number of atoms of the element in the molecule.
       If a mixture,such as concrete, pz could be the percent of the atoms
    """
    wa: float
    """Atomic weight"""


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
        self.calc_ge0_ge1()
        self.sigmas = {}

    def calc_ge0_ge1(self):
        """Define log intervals over the lower and upper cutoff energies"""
        self.ge1 = self.mge - 1 / log(self.up / self.ap)
        self.ge0 = 1 - self.ge1 * log(self.ap)

    def calc_sigmas(self, interaction, cross_sections):
        data = [0] * self.mge
        for element in self.elements:
            etmp, ftmp = cross_sections[element.z]

            # For pair or triplet, insert an extra data point for threshold energy
            # and change cross-sections (because ...?)
            if interaction in (Interaction.PAIR, Interaction.TRIPLET):
                eth = 2 * rm if interaction == Interaction.PAIR else 4 * rm
                ftmp = ftmp - 3 * np.log(1 - eth / np.exp(etmp))
                ftmp.insert(0, ftmp[0])
                etmp.insert(0, log(eth))

            for k in range(self.mge):
                gle = (k - 1 - self.ge0) / self.ge1
                exp_gle = exp(gle)
                # Check within bounds
                if not etmp[0] <= gle < etmp[-1]:  # not within
                    if interaction in (Interaction.PAIR, Interaction.TRIPLET):
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
                    kk = etmp.searchsorted(gle, side="right")
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

                if (
                    interaction in (Interaction.PAIR, Interaction.TRIPLET)
                    and exp_gle > eth
                ):
                    sig = sig * (1 - eth / exp_gle) ** 3

                data[k] += element.pz * sig

        self.sigmas[interaction] = data


if __name__ == "__main__":
    from egsnrc.hatch import DATA_DIR, get_xsection_table

    PHOTO = Interaction.PHOTOELECTRIC
    photo_data = get_xsection_table(DATA_DIR / "xcom_photo.data")
    # element_ta = Element(z=73, pz=1, wa=0)
    element_C = Element(z=6, pz=1, wa=0)
    medium = Medium("C", [element_C], rho=1, ap=0.01, up=50.0, mge=100)
    medium.calc_sigmas(PHOTO, photo_data)
    print(photo_data[6])
    print(medium.sigmas[PHOTO])
