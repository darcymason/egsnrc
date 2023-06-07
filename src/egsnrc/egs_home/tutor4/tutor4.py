from pathlib import Path
import os
from egsnrc import egsfortran
from egsnrc import watch
from egsnrc.electr import electr
from egsnrc.shower import shower
from egsnrc.util import for_E18, fort_hex  # for debugging vs mortran
import logging
import numpy  # cannot use `np` as is an EGS var!!
from math import log  # for calculate_tstep_...


# Get all common blocks
from egsnrc.commons import *

from egsnrc.calcfuncs import calc_tstep_from_demfp, compute_eloss, compute_eloss_g


init_done = False  # run init() only once, else crashes (unknown reason)


# Specific common blocks for this code
#  COMMON block geom --------
geom = egsfortran.geom
zbound = geom.zbound

#  COMMON block score --------
score = egsfortran.score
iwatch = score.iwatch
escore = score.escore


logger = logging.getLogger("egsnrc")  # XXX later `egsnrc`

HEN_HOUSE = Path(os.environ["HEN_HOUSE"])

# Define paths based on current location
# this python file should be in the user code directory under egs_home
#  - e.g. in tutor1, tutor2, etc.
HERE = Path(__file__).resolve().parent
EGS_HOME = HERE.parent
EGS_CONFIG = os.environ["EGS_CONFIG"]
USER_CODE = HERE.name
PEGS_FILE = "tutor_data"


def ausgab(iarg, **kwargs):
    """
    In this AUSGAB routine for tutor4, we score the energy deposited
    in the various regions. This amounts to the total energy
    reflected, deposited and transmitted by the slab.

    For IARG=0, an electron or photon step is about to occur and we
    score the energy deposited, if any. Note that only electrons
    deposit energy during a step, and due to our geometry, electrons
    only take steps in region 2 - however there is no need to check.
    For IARG=1,2 and 4, particles have been discarded for falling
    below various energy cutoffs and all their energy is deposited
    locally (in fact EDEP = particles kinetic energy).
    For IARG=3, we are discarding the particle since it is in
    region 1 or 3, so score its energy.
    """
    # logger.info("In Python ausgab")
    if iwatch > 0:
        # egsfortran.flush_output()
        watch.watch(iarg, iwatch, **kwargs)  # handles printouts of data
        # iwatch is passed in score

    if iarg <= 4:
        irl = ir[np - 1]  # pick up current region number  ** 0-based
        msg = " AUSGAB irl,edep"
        # logger.debug(f"{msg:<35}:{irl:3}{fort_hex(edep)}")
        escore[irl - 1] += edep


def print_info():
    print("egsfortran values", flush=True)
    print("-----------------", flush=True)
    for name in ("egs_home", "user_code", "pegs_file"):
        print(f"{name}: ", getattr(egsfortran.egs_io, name), flush=True)
    print("\nEnvironment", flush=True)
    print("-----------", flush=True)
    print(f"HEN_HOUSE={str(HEN_HOUSE)}", flush=True)
    print(f"EGS_HOME={str(EGS_HOME)}", flush=True)
    print(f"EGS_CONFIG={str(EGS_CONFIG)}", flush=True)


# ---------------------------------------------------------------------
# STEP 1:  USER-OVERRIDE-OF-EGSnrc-MACROS
# ---------------------------------------------------------------------
# REPLACE {$MXMED} WITH {1}  # only 1 medium in the problem(default 10)
# REPLACE {$MXREG} WITH {3}  # only 3 geometric regions (default 2000)
# REPLACE {$MXSTACK} WITH {15}"less than 15 particles on stack at once"

# # Define a common to pass information to the geometry routine HOWFAR
# REPLACE {;COMIN/GEOM/;} WITH {;COMMON/GEOM/ZBOUND;$REAL ZBOUND;}
# REPLACE {$CALL-HOWNEAR(#);} WITH {
#    ;{P1} = HOWNEAR({P1},X(NP),Y(NP),Z(NP),IRL);}
# # Define a COMMON for scoring in AUSGAB
# REPLACE {;COMIN/SCORE/;} WITH {
#    ;COMMON/SCORE/ESCORE(3),IWATCH; $INTEGER IWATCH; REAL*8 ESCORE;
# }


# ---------------------------------------------------------------------
# STEP 2 PRE-HATCH-CALL-INITIALIZATION
# ---------------------------------------------------------------------

# paths are done below so are set after egs_check_arguments


# --------------------------------------------------------------------
# egsfortran.egs_init()
def init(iwatch=1, high_prec=False):
    global init_done

    # import importlib
    # importlib.reload(egsfortran)
    if not init_done:
        # print("---Before setting pegs_file and user_code --", flush=True)
        # print_info()
        egsfortran.egs_set_defaults()
        egsfortran.egs_check_arguments()
        egsfortran.flush_output()

        # print("COMMON IO")
        # print("---------")
        # for name in dir(egsfortran.egs_io):
        #     if not name.startswith("_"):
        #         print(f'   {name} =', getattr(egsfortran.egs_io, name))

        egsfortran.egs_io.egs_home = f"{str(EGS_HOME) + '/':<128}"  # need trailing "/"
        egsfortran.egs_io.pegs_file = f"{PEGS_FILE:<256}"
        egsfortran.egs_io.user_code = f"{USER_CODE:<64}"

        print("\n---After setting pegs_file and user_code --", flush=True)
        print_info()

        egsfortran.egs_init1()
        # ----- end equiv of egs_init

    # Gotta be a better way, but for now this works.
    #  Blanking the third line because "NAI" is the default value in this array (??)
    # media is defined as $TYPE, with media[24,1].
    # $TYPE is macro'd to CHARACTER*4 for F77
    #
    # f2py message explains it:
    #  character*4 media(24,1) is considered as "character media(24,1,4)";
    # "intent(c)" is forced
    media[0, 0] = b"T   "
    media[1, 0] = b"A   "
    media[2, 0] = b"    "
    # print(media)

    # vacuum in regions 1 and 3, TA in region 2
    med[0] = med[2] = 0
    med[1] = 1

    # Note take 1 off indices for f2py 0-based in Python
    ecut[2 - 1] = 1.5  #    terminate electron histories at 1.5 MeV in the plate#
    pcut[2 - 1] = 0.1  #    terminate   photon histories at 0.1 MeV in the plate#
    #                only needed for region 2 since no transport elsewhere#
    #                ECUT is total energy = 0.989   MeV kinetic energy

    # ---------------------------------------------------------------------
    # STEP 3   HATCH-CALL
    # ---------------------------------------------------------------------

    if not init_done:
        logger.info("  Start tutor1\n\n CALL HATCH to get cross-section data\n")
        egsfortran.hatch()  #     pick up cross section data for TA
        #                data file must be assigned to unit 12

        egsfortran.flush_output()  # gfortran only - else doesn't print all lines

        logger.info(
            "\n knock-on electrons can be created and any electron followed down to\n"
            "                                       "
            f"{ae[0]-prm:8.3} MeV kinetic energy\n"
            " brem photons can be created and any photon followed down to      \n"
            "                                       "
            f"{ap[0]:8.3f} MeV"
            # Compton events can create electrons and photons below these cutoffs
        )  # OUTPUT AE(1)-PRM, AP(1);

    # ---------------------------------------------------------------------
    # STEP 4  INITIALIZATION-FOR-HOWFAR and HOWNEAR
    # ---------------------------------------------------------------------
    # Standard tutor4.mortran has zbound=0.1 (without D0, so single precision)
    # Here mimic that exact single precision number to get identical results
    # Even before end of history one, a min(tustep, tperp (from hownear)) is
    # different between Python and Fortran without this.
    point1_single_prec = float.fromhex("0x1.99999a0000000p-4")
    # OR, change tutor4.mortran to 0.1D0 and then compare with 0.1 below
    geom.zbound = 0.1  # point1_single_prec  #      plate is 1 mm thick

    # ---------------------------------------------------------------------
    # STEP 5  INITIALIZATION-FOR-AUSGAB
    # ---------------------------------------------------------------------
    # Print header for output - which is all AUSGAB does in this case
    # print("                 Kinetic Energy(MeV)  charge  angle w.r.t.Z axis-degrees")
    # for i in range(len(iausfl)):
    #     iausfl[i] = 1
    escore[:] = 0.0  # zero scoring array before starting

    watch.high_prec = high_prec  # if True, show more decimal places
    score.iwatch = iwatch  # This determines the type and amount of output
    # =1 => print info about each interaction
    # =2 => print info about same + each electron step
    # =4 => create a file to be displayed by EGS_Windows
    #  Note that these files can be huge
    # IWATCH 1 and 2 outputs to unit 6, 4 to unit 13

    egsfortran.flush_output()
    watch.watch(-99, iwatch)  # Initializes calls to AUSGAB for WATCH
    egsfortran.flush_output()  # if change above to egsfortran.watch(...)
    # ---------------------------------------------------------------------
    # STEP 6   DETERMINATION-OF-INICIDENT-PARTICLE-PARAMETERS
    # ---------------------------------------------------------------------
    # Define initial variables for 20 MeV beam of electrons incident
    # perpendicular to the slab

    init_done = True


def main(iqin=-1, iwatch=1, high_prec=False, ncase=10):
    # iqin here only to make generating validation data faster
    # The "in"s are local variables
    init(iwatch, high_prec)
    # et_control.exact_bca = False
    # et_control.spin_effects = True
    # iqin=-1  #                incident charge - electrons
    ein = 20.0 + prm
    ei = 20.0  #    20 MeV kinetic energy"
    xin = yin = zin = 0.0  #      incident at origin
    uin = vin = 0.0
    win = 1.0  #  moving along Z axis
    irin = 2  #                 starts in region 2, could be 1
    wtin = 1.0  #               weight = 1 since no variance reduction used
    # ---------------------------------------------------------------------
    # STEP 7   SHOWER-CALL
    # ---------------------------------------------------------------------
    # initiate the shower 10 times

    # INITIATE THE SHOWER NCASE TIMES
    #  ncase=10 # now from function call parameter

    callbacks = {"hownear": hownear, "howfar": howfar, "ausgab": ausgab}

    for i in range(ncase):
        if (iwatch != 0) and (iwatch != 4):
            logger.debug(
                "\n INITIAL SHOWER VALUES             :"
                f"    1{ei:9.3f}{iqin:4}{irin:4}"
                f"{xin:8.3f}{yin:8.3f}{zin:8.3f}"
                # below should be 8.3 like x,y,z but get extra spaces
                f"{uin:7.3f}{vin:7.3f}{win:7.3f}"
                f"{latchi:10}{wtin:10.3E}",
            )
        shower(iqin, ein, xin, yin, zin, uin, vin, win, irin, wtin, callbacks)

        # egsfortran.flush_output()
        watch.watch(-1, iwatch)  # print a message that this history is over
        # egsfortran.flush_output()

    # -----------------------------------------------------------------
    # STEP 8   OUTPUT-OF-RESULTS
    # -----------------------------------------------------------------

    # normalize to % of total input energy
    anorm = 100.0 / ((ein + iqin * prm) * ncase)
    total = sum(escore)

    msgs = (
        " Fraction of energy reflected from plate=",
        " Fraction of energy deposited in plate=",
        " Fraction of energy transmitted through plate=",
    )
    logger.info("\n")
    for i in range(3):
        logger.info(f"{msgs[i]:<49}{escore[i]*anorm:15.8f}%")

    logger.info(" " * 49 + "-" * 11)

    msg = " Total fraction of energy accounted for="
    logger.info(f"{msg:<49}{total*anorm:15.8f}%\n\n\n")

    # -----------------------------------------------------------------
    # STEP 9   finish run
    # -----------------------------------------------------------------
    egsfortran.egs_finish()

    # Expected hatch report for medium (line for pure Mortran/Fortran tutor1)
    # Medium            1
    # sige =    1.7946410123827996        1.7870816572288755
    # monotone =  T T

# *********************************************************************
#
def howfar():
    """Can particle can go distance ustep without crossing a boundary

    The following is a general specification of HOWFAR
    Given a particle at (X,Y,Z) in region IR and going in direction
    (U,V,W), this routine answers the question, can the particle go
    a distance USTEP without crossing a boundary?
            If yes, it merely returns
            If no, it sets USTEP=distance to boundary in the current
            direction and sets IRNEW to the region number   on the
            far side of the boundary (this can be messy in general!)

    The user can terminate a history by setting IDISC>0. Here we
    terminate all histories which enter region 3 or are going
    backwards in region 1

                    |               |
    REGION 1        |   REGION 2    |       REGION 3
                    |               |
    e- =========>   |               | e- or photon ====>
                    |               |
    vacuum          |     Ta        |       vacuum

    $REAL TVAL;
    COMIN/STACK,EPCONT,GEOM/;
           COMMON STACK contains X,Y,Z,U,V,W,IR and NP(stack pointer)
           COMMON EPCONT contains IRNEW, USTEP and IDISC
           COMMON GEOM contains ZBOUND
    """

    np_m1 = np - 1  # ** 0-based arrays

    if ir[np_m1] == 3:  # terminate this history: it is past the plate
        epcont.idisc = 1
        # logger.info("howfar  irl 3, idisc = 1")
        return

    if ir[np_m1] == 2:  # We are in the Ta plate - check the geometry
        if w[np_m1] > 0.0:  # going forward - consider first since  most frequent
            tval = (zbound - z[np_m1]) / w[np_m1]  # tval is dist to boundary
            #                                     in this direction
            # if tval > ustep, can take requested step, fall through to `return`
            if tval <= ustep:
                epcont.ustep = tval
                epcont.irnew = 3
        elif w[np_m1] < 0.0:  # going back towards origin
            tval = -z[np_m1] / w[np_m1]  # distance to plane at origin
            # if tval > ustep, can take requested step, fall through to `return`
            if tval <= ustep:
                epcont.ustep = tval
                epcont.irnew = 1
        # else w[np_m1] == 0.0, cannot hit boundary
        # logger.info("howfar region 2")
        return

    # Not region 3 or 2, must be 1, region with source
    if w[np_m1] > 0.0:  # this must be a source particle on z=0 boundary
        epcont.ustep = 0.0
        epcont.irnew = 2
        # logger.info("howfar region 1 going to 2")
    else:
        # it must be a reflected particle-discard it
        epcont.idisc = 1
        # logger.info("howfar reflected region 1, idisc=1")


def hownear(x, y, z, irl):
    """Determine distance to the closest boundary

    Given a particle at (x,y,z) in region irl, HOWNEAR answers the
    question, What is the distance to the closest boundary?

    In general this can be a complex subroutine.

    Parameters
    ----------
    x,y,z: REAL
        Position of particle in region irl
    irl: INTEGER
        Current region

    Returns
    -------
    REAL
        Distance to the closest boundary
    """
    # print(f"In Python hownear, with pos = ({x}, {y}, {z}), irl={irl}")

    if irl == 2:  # We are in the Ta plate - check the geometry
        return min(z, (zbound - z))  # 'terp'

    # else in region 1 or 3, can't return anything sensible
    raise ValueError(f"Called hownear in region {irl}")


if __name__ == "__main__":
    import sys
    import cProfile as profile
    import pstats

    HERE = Path(__file__).resolve().parent
    TEST_DATA = HERE.parent.parent / "tests" / "data"

    high_prec = False
    filename = HERE / "profile.stats"
    if len(sys.argv) == 1:
        # profile.run(
        #     "main(iwatch=iwatch, high_prec=high_prec, ncase=2000)",
        #     str(filename)
        # )
        # stats = pstats.Stats('profile.stats')
        main(iwatch=iwatch, high_prec=high_prec, ncase=2000)

        # Clean up filenames for the report
        # stats.strip_dirs()

        # Sort the statistics by the cumulative time spent
        # in the function
        # stats.sort_stats('cumulative')
        # stats.print_stats(20)
    elif sys.argv[1] == "gen":
        # Else, generating validation data for tests
        # generate for both e- and e+
        # filename = sys.argv[1]
        # for now, user still has to redirect to ../../tests/data/(filename)
        print("# e-   ----------------------------")
        main(-1, iwatch=iwatch, high_prec=high_prec)
        print("\n\n# e+   ----------------------------")
        main(+1, iwatch=iwatch, high_prec=high_prec)
    elif sys.argv[1] == "e+":
        main(iqin=+1, iwatch=2, high_prec=True, ncase=20)
    elif sys.argv[1] == "e-":
        main(iqin=-1, iwatch=2, high_prec=True, ncase=20)
    else:
        print(f"Unknown argument {sys.argv[1]}")
