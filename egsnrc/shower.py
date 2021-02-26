from egsnrc import egsfortran
from .commons import *
from .angles import uphi
from .electr import electr
import numpy

def shower(iqi, ei, xi, yi, zi, ui, vi, wi, iri, wti, callbacks):
    """Track a particle and others created during its simulation

    Parameters
    ----------
    iqi: int
        initial particle charge
    ei: float
        initial shower energy
    xi, yi, zi:  (float, float, float)
        initial co-ordinates
    ui, vi, wi: (float, float, float)
        initial direction cosines
    iri: int
        initial region number
    wti: float
        initial weight
    callbacks: dict
        Callback functions for `hownear`, `howfar` and `ausgab`

    """

    hownear = callbacks['hownear']
    howfar = callbacks['howfar']
    ausgab = callbacks['ausgab']
    egsfortran.ausgab = ausgab  # for functions still in egsfortran to call

    # msg = ", ".join(
    #     f"{x}={locals()[x]}"
    #     for x in "iqi,ei,xi,yi,zi,ui,vi,wi,iri,wti".split(",")
    # )
    # logger.info(f"Called shower with {msg}")
    ircode = numpy.array(0)  # meed rank-0 array for output var in f2py

    # $ comin_shower # DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
                    # COMIN/DEBUG,STACK,UPHIOT,RANDOM/
    # # Local variables
    # DOUBLE PRECISION
    # deg,    energy for pi-zero option
    # dpgl,   angle factor for pi-zero option
    # dei,    incident energy for pi-zero option
    # dpi,    intermediate factor for pi-zero option
    # dcsth,  random number for pi-zero option
    # dcosth, cos(theta) for pi-zero option
    # pi0msq  pi-zero mass squared (in MeV**2)

    # real*8
    # dneari, initial distance to closest boundary
    # csth,  random number for pi-zero option

    # integer*4
    # ircode  status returned by ELECTR or PHOTON

    pi0msq = 1.8215416  # PI-ZERO MASS (MEV) SQUARED

    stack.np=1
    stack.npold = stack.np # Set the old stack counter
    dneari=0.0
    iq[0]=iqi
    e[0]=ei
    u[0]=ui
    v[0]=vi
    w[0]=wi

    # TRANSFER PROPERTIES TO [0] FROM I  # ** 0-based, is to [1] in Mortran
    x[0]=xi
    y[0]=yi
    z[0]=zi
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
        # uphiot.costhe=dcosth
        # uphiot.sinthe=dsqrt(1.0-dcosth*dcosth)
        # iq[0]=0
        # e[0]=deg/2.
        # uphi(2,1)
        # stack.np=2
        # deg=dei-dpi*dcsth
        # dpgl=dpi-dei*dcsth
        # dcosth=dpgl/deg
        # uphiot.costhe=dcosth
        # uphiot.sinthe=-sqrt(1.0-dcosth*dcosth)
        # iq[2-1]=0
        # e[2-1]=deg/2.
        # uphi(3,2)


    while np > 0:
        #  DEFAULT FOR $ KERMA-INSERT; IS ; (NULL)
        if  iq[np-1] == 0:  # np-1 for ** 0-based in Python
            egsfortran.photon(ircode, howfar)
        else:
            # Note, callbacks have to be passed as extra parameters
            # even if not in the mortran call arguments,
            # unless intent(callback,hide) is used in f2py comments,
            # in which case, need to set `egsfortran.hownear = hownear`
            # egsfortran.electr(ircode, howfar) #, hownear,
            #     calc_tstep_from_demfp,
            #     compute_eloss, compute_eloss_g
            # )
            ircode = electr(hownear, howfar, ausgab)
        # egsfortran.flushoutput()
