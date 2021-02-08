
from pathlib import Path
import os
from egsnrc import egsfortran
import logging
import numpy  # cannot use `np` as is an EGS var!!
from math import log  # for calculate_tstep_...

# Get all common blocks
from egsnrc.commons import *


init_done = False  # run init() only once, else crashes (unknown reason)


# Specific common blocks for this code
#  COMMON block geom --------
geom = egsfortran.geom
zbound = geom.zbound

#  COMMON block score --------
score = egsfortran.score
iwatch = score.iwatch
escore = score.escore


logger = logging.getLogger('egsnrc')  # XXX later `egsnrc`

HEN_HOUSE = Path(os.environ["HEN_HOUSE"])
EGS_HOME = Path(os.environ['EGS_HOME'])


# ******************************************************************
#                                National Research Council of Canada
def shower(iqi,ei,xi,yi,zi,ui,vi,wi,iri,wti):
    #
    # ******************************************************************

    # stack
    global e, x, y, z, u, v, w, dnear, wt, iq, ir, latch, latchi, np, npold
    # uphiot
    global theta, sinthe, costhe, sinphi, cosphi, pi, twopi, pi5d2

    # msg = ", ".join(f"{x}={locals()[x]}" for x in "iqi,ei,xi,yi,zi,ui,vi,wi,iri,wti".split(","))
    # logger.info(f"Called shower with {msg}")
    ircode = numpy.array(0)  # meed rank-0 array for output var in f2py

    # $ comin_shower # DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
                    # COMIN/DEBUG,STACK,UPHIOT,RANDOM/

    # # Input variables
    # ;real*8 EI,      # initial shower energy
    #       XI,YI,ZI,# initial co-ordinates
    #       UI,VI,WI,# initial direction cosines
    #       WTI # initial weight

    # ;integer*4
    #       IQI,     # initial particle charge
    #       IRI # initial region number

    # # Local variables
    # DOUBLE PRECISION
    #       DEG,    # energy for pi-zero option
    #       DPGL,   # angle factor for pi-zero option
    #       DEI,    # incident energy for pi-zero option
    #       DPI,    # intermediate factor for pi-zero option
    #       DCSTH,  # random number for pi-zero option
    #       DCOSTH, # cos(theta) for pi-zero option
    #       PI0MSQ # pi-zero mass squared (in MeV**2)

    # ;real*8 DNEARI, # initial distance to closest boundary
    #       CSTH # random number for pi-zero option

    # ;integer*4
    #       IRCODE # status returned by ELECTR or PHOTON

    PI0MSQ = 1.8215416  # PI-ZERO MASS (MEV) SQUARED

    stack.np=1
    stack.npold = stack.np # Set the old stack counter
    dneari=0.0
    iq[0]=iqi
    e[0]=ei
    u[0]=ui; v[0]=vi; w[0]=wi


    # TRANSFER PROPERTIES TO [0] FROM I
    x[0]=xi; y[0]=yi; z[0]=zi
    ir[0]=iri
    wt[0]=wti
    dnear[0]=dneari
    latch[0]=latchi

    if iqi == 2:
        # PI-ZERO OPTION
        # if EI <= PI0MSQ) [OUTPUT EI;    corrected Oct 24 1995 e-mail Hideo H
        #                   noted by      Dr.  Muroyama at Nagoya University
        raise NotImplementedError("egsnrc Python not tested for PI-ZERO")
        # if ei**2 <= pi0msq:
        #     msg = (
        #         ' Stopped in subroutine SHOWER---PI-ZERO option invoked'
        #         f' but the total energy was too small (EI={ei} MeV)'
        #     )
        #     raise ValueError(msg)

        # csth = randomset()
        # dcsth=csth; dei=ei; dpi=dsqrt(dei*dei-pi0msq)
        # deg=dei+dpi*dcsth; dpgl=dpi+dei*dcsth; dcosth=dpgl/deg
        # costhe=dcosth; sinthe=dsqrt(1.0-dcosth*dcosth)  # "1.D0" -> "1.0" ???
        # iq[0]=0; e[0]=deg/2.
        # egsfortran.uphi(2,1)
        # stack.np=2
        # deg=dei-dpi*dcsth; dpgl=dpi-dei*dcsth; dcosth=dpgl/deg
        # costhe=dcosth; sinthe=-dsqrt(1.0-dcosth*dcosth)  # "1.D0" -> "1.0" ???
        # iq[2-1]=0; e[2-1]=deg/2.
        # egsfortran.uphi(3,2)


    while np > 0:
        #  DEFAULT FOR $ KERMA-INSERT; IS ; (NULL)
        if  iq[np-1] == 0:  # -1 for ** 0-based in Python
            egsfortran.photon(ircode, howfar)
        else:
            # Note, callbacks have to be passed as extra parameters
            # even if not in the mortran call arguments,
            # unless intent(callback,hide) is used in f2py comments,
            # in which case, need to set `egsfortran.hownear = hownear`
            egsfortran.electr(ircode, howfar, hownear,
                calc_tstep_from_demfp,
                compute_eloss, compute_eloss_g
            )
        # egsfortran.flushoutput()
    # ---------------- end of subroutine shower


# Define paths based on current location
# this python file should be in the user code directory under egs_home
#  - e.g. in tutor1, tutor2, etc.
HERE  = Path(__file__).resolve().parent
EGS_HOME = HERE.parent
EGS_CONFIG = os.environ['EGS_CONFIG']
USER_CODE = HERE.name
PEGS_FILE = "tutor_data"


def print_info():
    print("egsfortran values", flush=True)
    print("-----------------", flush=True)
    for name in ('egs_home', 'user_code', 'pegs_file'):
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
def init():
    global init_done

    if init_done:
        return
    # print("---Before setting pegs_file and user_code --", flush=True)
    # print_info()

    egsfortran.egs_set_defaults()
    egsfortran.egs_check_arguments()
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
    # media is defined as $TYPE, with media[24,1]. $TYPE is macro'd to CHARACTER*4 for F77
    #
    # f2py message explains it:
    #  character*4 media(24,1) is considered as "character media(24,1,4)"; "intent(c)" is forced
    media[0,0] = b'T   '
    media[1,0] = b'A   '
    media[2,0] = b'    '
    # print(media)

    # vacuum in regions 1 and 3, TA in region 2
    med[0] = med[2] = 0
    med[1] = 1

    # Note take 1 off indices for f2py 0-based in Python
    ecut[2-1]=1.5;  #    terminate electron histories at 1.5 MeV in the plate#
    pcut[2-1]=0.1;  #    terminate   photon histories at 0.1 MeV in the plate#
    #                only needed for region 2 since no transport elsewhere#
    #                ECUT is total energy = 0.989   MeV kinetic energy


    # ---------------------------------------------------------------------
    # STEP 3   HATCH-CALL
    # ---------------------------------------------------------------------

    logger.info('  Start tutor1\n\n CALL HATCH to get cross-section data\n')
    egsfortran.hatch()  #     pick up cross section data for TA
    #                data file must be assigned to unit 12

    # egsfortran.flushoutput()  # gfortran only - else doesn't print all lines

    logger.info(
        ' knock-on electrons can be created and any electron followed down to\n'
        "                                        "
        f'{ae[0]-prm:8.3} MeV kinetic energy\n'
        ' brem photons can be created and any photon followed down to      \n'
        "                                        "
        f'{ap[0]:8.3} MeV '
        # Compton events can create electrons and photons below these cutoffs
    )# OUTPUT AE(1)-PRM, AP(1);

    # ---------------------------------------------------------------------
    # STEP 4  INITIALIZATION-FOR-HOWFAR and HOWNEAR
    # ---------------------------------------------------------------------
    geom.zbound=0.1  #      plate is 1 mm thick

    # ---------------------------------------------------------------------
    # STEP 5  INITIALIZATION-FOR-AUSGAB
    # ---------------------------------------------------------------------
    # Print header for output - which is all AUSGAB does in this case
    # print("                 Kinetic Energy(MeV)  charge  angle w.r.t.Z axis-degrees")
    for i in range(3):
        escore[i] = 0.0  # zero scoring array before starting

    score.iwatch=1  # This determines the type and amount of output
                    # =1 => print info about each interaction
                    # =2 => print info about same + each electron step
                    # =4 => create a file to be displayed by EGS_Windows
                    #  Note that these files can be huge
                    # IWATCH 1 and 2 outputs to unit 6, 4 to unit 13

    egsfortran.watch(-99,iwatch);   # Initializes calls to AUSGAB for WATCH
    # ---------------------------------------------------------------------
    # STEP 6   DETERMINATION-OF-INICIDENT-PARTICLE-PARAMETERS
    # ---------------------------------------------------------------------
    # Define initial variables for 20 MeV beam of electrons incident
    # perpendicular to the slab

    init_done = True


def main(iqin=-1):  # iqin here only to make generating validation data faster
    # The "in"s are local variables
    init()
    # et_control.exact_bca = False
    # et_control.spin_effects = True
    # iqin=-1  #                incident charge - electrons
    ein=20 + prm
    ei=20.0  #    20 MeV kinetic energy"
    xin = yin = zin = 0.0  #      incident at origin
    uin = vin = 0.0; win=1.0  #  moving along Z axis
    irin=2  #                 starts in region 2, could be 1
    wtin=1.0  #               weight = 1 since no variance reduction used
    # ---------------------------------------------------------------------
    # STEP 7   SHOWER-CALL
    # ---------------------------------------------------------------------
    # initiate the shower 10 times

    ncase=10  # INITIATE THE SHOWER NCASE TIMES

    for i in range(ncase):
        if (iwatch != 0) and (iwatch != 4):
            print(
            "\n INITIAL SHOWER VALUES             :"
            f"    1{ei:9.3f}{iqin:4}{irin:4}"
            f"{xin:8.3f}{yin:8.3f}{zin:8.3f}"
            f"{uin:7.3f}{vin:7.3f}{win:7.3f}"  # should be 8.3 like x,y,z but get extra spaces
            f"{latchi:10}{wtin:10.3E}"
            )
            shower(iqin,ein,xin,yin,zin,uin,vin,win,irin,wtin)
            egsfortran.watch(-1,iwatch)  # print a message that this history is over

    # -----------------------------------------------------------------
    # STEP 8   OUTPUT-OF-RESULTS
    # -----------------------------------------------------------------
    # note output is at the end of each history in subroutine ausgab

    # -----------------------------------------------------------------
    # STEP 9   finish run
    # -----------------------------------------------------------------
    egsfortran.egs_finish()

    # Expected hatch report for medium (line for pure Mortran/Fortran tutor1)
    # Medium            1  sige =    1.7946410123827996        1.7870816572288755       monotone =  T T

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
        return

    # Not region 3 or 2, must be 1, region with source
    if w[np_m1] >  0.0:  # this must be a source particle on z=0 boundary
        epcont.ustep = 0.0
        epcont.irnew = 2
    else:
        # it must be a reflected particle-discard it
        epcont.idisc = 1


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
        return min(z, (zbound - z) )  # 'terp'

    # else in region 1 or 3, can't return anything sensible
    raise ValueError(f'Called hownear in region {irl}')


def compute_drange(lelec, medium, eke1, eke2, lelke1, elke1, elke2):
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
    fedep = 1 - eke2/eke1

    # evaluate the logarithm of the midpoint energy
    elktmp = 0.5*(elke1+elke2+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)))

    # *** -1 for 0-based in Python
    lelktmp = lelke1 - 1 # was = lelke1
    medium -= 1

    if lelec < 0:
        # $EVALUATE dedxmid USING ededx(elktmp)
        dedxmid = ededx1[lelktmp,medium]*elktmp+ ededx0[lelktmp,medium]
        dedxmid = 1/dedxmid
        aux = ededx1[lelktmp,medium]*dedxmid
        #  aux = ededx1(lelktmp,medium)/dedxmid"
    else:
        # $EVALUATE dedxmid USING pdedx(elktmp)
        dedxmid = pdedx1[lelktmp,medium]*elktmp+ pdedx0[lelktmp,medium]
        dedxmid = 1/dedxmid
        aux = pdedx1[lelktmp,medium]*dedxmid
        #  aux = pdedx1(lelktmp,medium)/dedxmid

    aux = aux*(1+2*aux)*(fedep/(2-fedep))**2/6

    return fedep*eke1*dedxmid*(1+aux)


def calc_tstep_from_demfp(qel,lelec, medium, lelke, demfp, sig, eke, elke, total_de):
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

    # print("fn:", ",".join(str(x) for x in (qel,lelec, medium, demfp, sig, eke, elke, total_de)) )
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


def compute_eloss(lelec, medium, step, eke, elke, lelke):
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
    de = dedxmid*step*rhof # IK: rhof scaling bug, June 9 2006
                            # rhof scaling must be done here and NOT in
                            # $ COMPUTE-ELOSS-G
    fedep = de / eke
    de = de*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-0.25*fedep*(2-aux*(4-aux)))))

    return de


def compute_eloss_g(lelec, medium, step, eke, elke, lelke, range_):
    """A generalized version of `compute_eloss`"""
    # ** 0-based arrays in Python
    medium_m1 = medium - 1
    lelke_m1 = lelke - 1

    # Note: range_ep IS 0-based already in first dimn

    qel = 0 if lelec==-1 else 1  # recalc here to not bother passing in both
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
            de = eke - te[medium_m1]*0.99
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

    return de


def calculate_xi(lelec, medium, ekems, rmt2, rmsq, xccl, blccl, step):
    # ** 0-based arrays
    medium_m1 = medium - 1

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


if __name__ == "__main__":
    import sys
    HERE  = Path(__file__).resolve().parent
    TEST_DATA = HERE.parent.parent / "tests" / "data"

    if len(sys.argv) == 1:
        main()
    else:
        # Else, generating validation data for tests
        # generate for both e- and e+
        # filename = sys.argv[1]
        # for now, user still has to redirect to ../../tests/data/(filename)
        if sys.argv[1] != "gen":
            print("Only accept 'gen' as optional argument")
        else:
            print("# e-   ----------------------------")
            main(-1)
            print("\n\n# e+   ----------------------------")
            main(+1)

