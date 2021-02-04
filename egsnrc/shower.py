from egsnrc2py import egsfortran

ranlux = egsfortran.ranlux
rng_seed = egsfortran.randomm.rng_seed
rng_array = egsfortran.randomm.rng_array


np = egsfortran.stack.np
npold = egsfortran.stack.npold

# CALLBACKS ---- 

def randomset():
# randomm
    global rng_seed, rng_array, seeds

    if rng_seed > 24:
        ranlux(rng_array)
        rng_seed = 1

    random_num = rng_array[rng_seed-1]
    rng_seed += 1

    return random_num


# ******************************************************************
#                                National Research Council of Canada
def SHOWER(iqi,ei,xi,yi,zi,ui,vi,wi,iri,wti):
    #                                                                   
    # ******************************************************************

    # stack
    global e, x, y, z, u, v, w, dnear, wt, iq, ir, latch, latchi, np, npold
    # uphiot
    global theta, sinthe, costhe, sinphi, cosphi, pi, twopi, pi5d2

    
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

    PI0MSQ = 1.8215416D4  # PI-ZERO MASS (MEV) SQUARED

    np=1
    npold = np # Set the old stack counter
    dneari=0.0
    iq[0]=iqi; e[0]=ei; u[0]=ui; v[0]=vi; w[0]=wi


    # TRANSFER PROPERTIES TO [0] FROM I
    x[0]=xi;y[0]=yi;z[0]=zi;ir[0]=iri
    wt[0]=wti;dnear[0]=dneari;latch[0]=latchi

    if IQI == 2:
        # PI-ZERO OPTION
        # IF(EI <= PI0MSQ) [OUTPUT EI;    corrected Oct 24 1995 e-mail Hideo H 
        #                   noted by      Dr.  Muroyama at Nagoya University
        if ei**2 <= pi0msq:
            msg = (
                ' Stopped in subroutine SHOWER---PI-ZERO option invoked'
                f' but the total energy was too small (EI={ei} MeV)'
            )
            raise ValueError(msg)

        csth = randomset()
        dcsth=csth; dei=ei; dpi=dsqrt(dei*dei-pi0msq)
        deg=dei+dpi*dcsth; dpgl=dpi+dei*dcsth; dcosth=dpgl/deg
        costhe=dcosth; sinthe=dsqrt(1.d0-dcosth*dcosth)
        iq[0]=0; e[0]=deg/2.
        call uphi(2,1)
        np=2
        deg=dei-dpi*dcsth; dpgl=dpi-dei*dcsth; dcosth=dpgl/deg
        costhe=dcosth; sinthe=-dsqrt(1.d0-dcosth*dcosth)
        iq(2)=0; e(2)=deg/2.
        egsfortran.uphi(3,2)


    while np > 0:
        #  DEFAULT FOR $ KERMA-INSERT; IS ; (NULL) 
        if  iq[np-1] == 0:
            call egsfortran.photon(ircode)
        else:
            call egsfortran.electr(ircode)

    # end of subroutine shower

