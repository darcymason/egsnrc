# media.py
"""Classes for Media and Elements"""
from dataclasses import dataclass

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


@dataclass  # XXX make frozen=True later?
class Medium:
    """Contain information for a specific physical medium"""
    name: str
    elements: list[Element]
    rho: float

