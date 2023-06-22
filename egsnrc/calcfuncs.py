from math import log
from egsnrc.commons import *
from egsnrc import config

import logging
logger = logging.getLogger("egsnrc")

def compute_drange(lelec: int, medium: int, eke1: float, eke2: float, lelke1: int, elke1: float, elke2:float) -> float:
    """Computes path-length traveled going from energy `eke1` to `eke2`

    both energies being in the same interpolation bin,
    given by `lelke1`. `elke1` and `elke2` are the logarithms of
    'eke1' and `eke2`. The expression is based on logarithmic interpolation as
    used in EGSnrc (i.e. dedx = a + b*Log(E) ) and a power series expansion
    of the ExpIntegralEi function that is the result of the integration.

    Parameters
    ----------
    lelec: INTEGER
        Charge of the particle, either -1 or +1

    medium: INTEGER
        Current medium

    Returns
    -------
    REAL
        path-length traveled going from energy `eke1` to `eke2`

    """
    fedep: float = 1 - eke2/eke1

    # evaluate the logarithm of the midpoint energy
    elktmp: float = 0.5*(elke1+elke2+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)))

    # *** -1 for 0-based in Python
    lelktmp = lelke1 - 1 # was = lelke1
    medium_m1 = medium - 1

    if lelec < 0:
        # $EVALUATE dedxmid USING ededx(elktmp)
        dedxmid = ededx1[lelktmp,medium_m1]*elktmp+ ededx0[lelktmp,medium_m1]
        dedxmid = 1/dedxmid
        aux = ededx1[lelktmp,medium_m1]*dedxmid
        #  aux = ededx1(lelktmp,medium_m1)/dedxmid"
    else:
        # $EVALUATE dedxmid USING pdedx(elktmp)
        dedxmid = pdedx1[lelktmp,medium_m1]*elktmp+ pdedx0[lelktmp,medium_m1]
        dedxmid = 1/dedxmid
        aux = pdedx1[lelktmp,medium_m1]*dedxmid
        #  aux = pdedx1(lelktmp,medium_m1)/dedxmid

    aux = aux*(1+2*aux)*(fedep/(2-fedep))**2/6

    return fedep*eke1*dedxmid*(1+aux)


def calc_tstep_from_demfp(qel,lelec, medium: int, lelke, demfp, sig, eke, elke, total_de):
    """Calculate path length to the next discrete interaction

    Once the sub-threshold processes energy loss to the next discrete
    interaction is determined, the corresponding path-length has to be
    calculated. This is done by this function. This function
    assumes the energy at the begining to be `eke`, the logarithm of it
    `elke`, `lelke` - the corresponding interpolation index and makes
    use of `compute_drange`.
    """
    # in: medium, qel, leklef (for compute-drange)
    # in/out:  compute_tstep, total_tstep,
    # global e_array, epcont.eke, epcont.elke, bounds.vacdst, eke0[], eke1[]

    # print("fn:", ",".join(str(x) for x in (qel,lelec: int, medium, demfp, sig, eke, elke, total_de)) )
    fedep = total_de
    ekef  = eke - fedep

    # *** 0-based array
    medium_m1 = medium - 1
    if  ekef <= e_array[1-1,medium_m1]:
        tstep = vacdst
    else:
        elkef = log(ekef)
        # Unhandled macro '$ SET INTERVAL elkef,eke;'->
        # Below line from tutor4_linux.f
        lelkef=int(eke1[medium_m1]*elkef+eke0[medium_m1])  # XXX note fortran had implicit real->int conversion
        if  lelkef == lelke:
            #  initial and final energy are in the same interpolation bin
            # --- Inline replace: $ COMPUTE_DRANGE(eke,ekef,lelke,elke,elkef,tstep); -----
            tstep = compute_drange(lelec, medium, eke,ekef,lelke,elke,elkef)
        else:
            #  initial and final energy are in different interpolation bins,
            #  calc range from ekef to E(lelkef+1) and from E(lelke) to eke
            #  and add the pre-calculated range from E(lelkef+1) to E(lelke)
            ekei = e_array[lelke-1,medium_m1]
            elkei = (lelke - eke0[medium_m1])/eke1[medium_m1]
            # --- Inline replace: $ COMPUTE_DRANGE(eke,ekei,lelke,elke,elkei,tuss); -----
            tuss = compute_drange(lelec, medium, eke,ekei,lelke,elke,elkei)
            ekei = e_array[lelkef+1-1,medium_m1]  # 0-based -1
            elkei = (lelkef + 1 - eke0[medium_m1])/eke1[medium_m1]
            # --- Inline replace: $ COMPUTE_DRANGE(ekei,ekef,lelkef,elkei,elkef,tstep); -----
            tstep = compute_drange(lelec, medium, ekei,ekef,lelkef,elkei,elkef)
            # Note: range_ep IS 0-based already in first dimn
            tstep=tstep+tuss+range_ep[qel,lelke-1,medium_m1]-range_ep[qel,lelkef+1-1,medium_m1]

    return tstep



# Get exact match to Fortan constants (single prec) used in EGSnrc Mortran
# Python is always double prec so it represented this constant differently
# than Fortan without the "D0" suffix to mark it as double
# if `test_precision`, then use the EGSnrc Mortran to match exactly for
# comparing outputs exactly
if config.test_precision:
    point_333333 = float.fromhex('0x1.55553e0000000p-2')
    point_99 =float.fromhex('0x1.fae1480000000p-1')
else:
    point_333333 = 0.333333
    point_99 = 0.99


def compute_eloss(lelec: int, medium: int, step, eke, elke, lelke):
    """"Compute the energy loss due to sub-threshold processes for a path-length `step`.

    The energy at the beginning of the step is `eke`, `elke`=log(`eke`),
    `lelke` is the interpolation index.
    The formulae are based on the logarithmic interpolation for dedx
    used in EGSnrc.

    Returns
    -------
    REAL
        energy loss

    Note
    ----
    Assumes that initial and final energy are in the same interpolation bin.

    """
    # print("fn: ",lelec, medium, step, eke, elke, lelke)

    # logger.debug(f"in compute-eloss:{fort_hex([step, eke, elke])}{lelke:4}")
    # ** 0-based
    medium_m1 = medium - 1
    lelke_m1 = lelke - 1

    if lelec < 0:
        dedxmid = ededx1[lelke_m1, medium_m1]*elke+ ededx0[lelke_m1, medium_m1]  # EVALUATE dedxmid USING ededx(elke)
        aux = ededx1[lelke_m1, medium_m1]/dedxmid
    else:
        dedxmid = pdedx1[lelke_m1, medium_m1]*elke+ pdedx0[lelke_m1, medium_m1]  # EVALUATE dedxmid USING pdedx(elke)
        aux = pdedx1[lelke_m1, medium_m1]/dedxmid

    # de = dedxmid*tuss #  Energy loss using stopping power at the beginning
    de = dedxmid*step*rhof  # IK: rhof scaling bug, June 9 2006
                            # rhof scaling must be done here and NOT in
                            # $ COMPUTE-ELOSS-G
    fedep = de / eke

    de *= 1 - 0.5*fedep*aux*(
        1 - point_333333*fedep*(aux - 1 -0.25*fedep*(2-aux*(4-aux)))
    )

    # logger.debug(f"out compute-eloss:{fort_hex(de)}")
    return de


def compute_eloss_g(lelec: int, medium: int, step: float, eke:float, elke:float, lelke:int, range_:float):
    """A generalized version of `compute_eloss`"""
    # ** 0-based arrays in Python
    medium_m1: int = medium - 1
    lelke_m1: int = lelke - 1

    # logger.debug(f"in compute-eloss-g:{fort_hex([step, eke, elke])}{lelke:4}")
    # Note: range_ep IS 0-based already in first dimn

    qel: int = 0 if lelec==-1 else 1  # recalc here to not bother passing in both
    tuss = range_ - range_ep[qel,lelke_m1,medium_m1] / rhof
        #  here tuss is the range between the initial energy and the next lower
        #  energy on the interpolation grid
    if tuss >= step:
        #  Final energy is in the same interpolation bin
        # --- Inline replace: $ COMPUTE_ELOSS(tustep,eke,elke,lelke,de); -----
        de = compute_eloss(lelec, medium, step, eke, elke, lelke)

    else:  # Must find first the table index where the step ends using
           #  pre-calculated ranges
        lelktmp = lelke
        tuss = (range_ - step) * rhof
        #  now tuss is the range of the final energy electron
        #  scaled to the default mass density from PEGS4

        if tuss <= 0:
            de = eke - te[medium_m1]*point_99
            #  i.e., if the step we intend to take is longer than the particle
            #  range, the particle energy goes down to the threshold
            # (eke is the initial particle energy)
            # originally the entire energy was lost, but msdist_xxx is not prepared
            # to deal with such large eloss fractions => changed July 2005.
        else:
            while tuss < range_ep[qel, lelktmp-1, medium_m1]:
                lelktmp -= 1
            lelktmp_m1 = lelktmp - 1  # *** 0-based arrays in Python
            elktmp = (lelktmp + 1 - eke0[medium_m1]) / eke1[medium_m1]
            eketmp = e_array[lelktmp_m1+1, medium_m1]
            # tuss = range_ep(qel,lelktmp+1,medium_m1) - tuss
            # IK: rhof scaling bug, June 9 2006: because of the change in
            #     compute_eloss above, we must scale tuss by rhof
            tuss = (range_ep[qel, lelktmp_m1+1, medium_m1] - tuss) / rhof
            # --- Inline replace: $ COMPUTE_ELOSS(tuss,eketmp,elktmp,lelktmp,de); -----
            de = compute_eloss(lelec, medium, tuss, eketmp, elktmp, lelktmp)
            de = de + eke - eketmp
    # logger.debug(f"out compute-eloss-g:{fort_hex(de)}")
    return de


def calculate_xi(lelec: int, medium: int, ekems: float, rmt2:float, rmsq:float, xccl:float, blccl:float, step:float):
    # ** 0-based arrays
    medium_m1: int = medium - 1

    p2 = ekems*(ekems+rmt2)
    beta2 = p2/(p2 + rmsq)
    chia2 = xccl/(4*blccl*p2)
    # Note that our chia2 is Moliere chia2/4
    # Note also that xcc is now old egs xcc**2
    xi = 0.5*xccl/p2/beta2*step
    if spin_effects:
        elkems = log(ekems)
        lelkems=int(eke1[medium_m1]*elkems+eke0[medium_m1])  # $ SET INTERVAL elkems,eke
        lelkems_m1 = lelkems - 1  # ** 0-based
        if lelec < 0:
            # EVALUATE etap USING etae_ms(elkems)
            etap = etae_ms1[lelkems_m1, medium_m1]*elkems+ etae_ms0[lelkems_m1, medium_m1]
            # EVALUATE xi_corr USING q1ce_ms(elkems)
            xi_corr = q1ce_ms1[lelkems_m1, medium_m1]*elkems+ q1ce_ms0[lelkems_m1, medium_m1]
        else:
            # EVALUATE etap USING etap_ms(elkems)
            etap = etap_ms1[lelkems_m1, medium_m1]*elkems+ etap_ms0[lelkems_m1, medium_m1]
            # EVALUATE xi_corr USING q1cp_ms(elkems)
            xi_corr = q1cp_ms1[lelkems_m1, medium_m1]*elkems+ q1cp_ms0[lelkems_m1, medium_m1]

        chia2 = chia2*etap
        xi = xi*xi_corr
        # EVALUATE ms_corr USING blcce(elkems)
        ms_corr = blcce1[lelkems_m1, medium_m1]*elkems+ blcce0[lelkems_m1, medium_m1]
        blccl = blccl*ms_corr  # not used after in this fn.  Needs to be returned?
    else:
        xi_corr = 1
        etap = 1

    xi = xi*(log(1+1./chia2)-1/(1+chia2))

    return xi, blccl
