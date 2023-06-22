# angles.py
from math import sqrt, sin, cos
from egsnrc.egsrandom import random_kfloat
from egsnrc.constants import TWO_PI
from egsnrc.params import SMALL_POLAR_ANGLE_THRESHOLD
from egsnrc import config

@config.device_jit
def uphi(rng_states, gid, p, sinthe, costhe):
    """Uniform phi distribution

    Returns
    -------
    u, v, w direction cosines
    """
    # In original Mortran UPHI(ientry, lvl):
    # SAVEd between calls: CTHET,PHI,CPHI,A,B,C,SINPS2,SINPSI,US,VS,SINDEL,COSDEL
    # ientry:
    #   1: calculate sinthe, costhe from theta (falls through to 2)
    #        (only used once in code, in PAIR after SET-PAIR-ANGLE, which
    #         itself calls uphi(2,1) for 1st particle, uphi(3,2) on -sinthe for 2nd),
    #         ?? ** but does not set theta itself??
    #      Fortran code for this step:
    #          SINTHE = sin(THETA)
    #          CTHET  = PI5D2-THETA  # pi5d2 = 5/2*PI, i.e. 2.5*PI
    #          COSTHE = sin(CTHET)
    #   2: Randomly set phi (falls through to 3)
    #   3: Adjust particle directions accordingly:
    #     Lvl 
    #         2: new particle: copy n-1 particle, use saved abc, ADJUST
    #         
    #         In old code source, only lvl 2 was used, i.e. (3,2).
    #         For reference, others were:
    #         1: old particle: abc = uvw (top particle), then ADJUST
    #         3: brems-gamma: abc=electron(n-1) uvw, ADJUST nth particle
    # :ADJUST: 
    #   calc uvw from theta, phi: 
    #       "a rotation is performed to transform direction cosines of the
    #       particle back to the physical frame (from the transport frame)"

    # unlike original Mortran calculate phi directly.
    # should be fast enough on todays computers
    phi = TWO_PI * random_kfloat(rng_states, gid)

    sinps2 = p.u * p.u + p.v * p.v
    # Small polar angle case
    if sinps2 < SMALL_POLAR_ANGLE_THRESHOLD:
        u = sinthe * cos(phi)
        v = sinthe * sin(phi)
        w = p.w * costhe
    else:
        sinpsi = sqrt(sinps2)
        us = sinthe * cos(phi)
        vs = sinthe * sin(phi)
        sindel = p.v / sinpsi
        cosdel = p.u / sinpsi

        u = p.w * cosdel * us - sindel * vs + p.u * costhe
        v = p.w * sindel * us + cosdel * vs + p.v * costhe
        w = -sinpsi * us + p.w * costhe
    return u, v, w



def UPHI(IENTRY,LVL):
    # "                                                                  "
    # "******************************************************************"
    # "   UPHI STANDS FOR 'UNIFORM PHI DISTRIBUTION'.                    "
    # "   SET COORDINATES FOR NEW PARTICLE OR RESET DIRECTION COSINES OF "
    # "   OLD ONE.  GENERATE RANDOM AZIMUTH SELECTION AND REPLACE THE    "
    # "   DIRECTION COSINES WITH THEIR NEW VALUES.                       "
    "******************************************************************"


    "Input variables"
    integer  IENTRY,LVL; "entry switches"

    "Local variables"
    $REAL CTHET,  "5/2*pi-THETA, used to evaluate cos(THETA) using the sine table"
        RNNO38, "random number for azimuthal angle selection"
        PHI,    "azimuthal scattering angle"
        CPHI,   "5/2*pi-PHI"
        A,B,C,  "direction cosines before rotation"
        SINPS2, "SINPS2=A*A+B*B"
        SINPSI, "Sqrt(SINPS2)"
        US,VS,  "x- and y- component of scattering vector"
        SINDEL,COSDEL;
                "aux. variables for the rotation algorithm"

    $INTEGER
        IARG,   "index for AUSGAB"
        LPHI,LTHETA,LCTHET,LCPHI;
                "indeces for sine table"

    $DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;
    save CTHET,PHI,CPHI,A,B,C,SINPS2,SINPSI,US,VS,SINDEL,COSDEL;

    $AUSCALL($UPHIAUSB);
    GO TO (:UPHI:,:UPHI2:,:NRK:),IENTRY;
    "   IENTRY OUT-OF-BOUNDS IF HERE"  GO TO :ERROR:;

    :UPHI:; "NOTE: AFB 88/12/12 ADDED SEMI-COLON, ELSE BUG WHEN OVERRIDING SIN"
            "TABLE LOOK-UP"
    $SET INTERVAL THETA,SINC;
    $EVALUATE SINTHE USING SIN(THETA);
    CTHET=PI5D2-THETA;$SET INTERVAL CTHET,SINC;
    $EVALUATE COSTHE USING SIN(CTHET);

    "   USE THE FOLLOWING ENTRY IF SINTHE AND COSTHE ARE ALREADY KNOWN."
    "   SELECT PHI UNIFORMLY OVER THE INTERVAL (0,TWO PI). THEN USE    "
    "   PWLF OF SIN FUNCTION TO GET SIN(PHI) AND COS(PHI).  THE COSINE "
    "   IS GOTTEN BY COS(PHI)=SIN(9*PI/4 - PHI).                       "

    :UPHI2:;

    " It is much faster to use the box method for azimuthal angle selection"
    " than the following                                                   "
    " $RANDOMSET RNNO38;
    " PHI=RNNO38*TWOPI;$SET INTERVAL PHI,SINC;
    " $EVALUATE SINPHI USING SIN(PHI);
    " CPHI=PI5D2-PHI;$SET INTERVAL CPHI,SINC;
    " $EVALUATE COSPHI USING SIN(CPHI);
    $SELECT-AZIMUTHAL-ANGLE(cosphi,sinphi);

    "   USE THE FOLLOWING ENTRY FOR THE SECOND OF TWO PARTICLES WHEN WE"
    "   KNOW TWO PARTICLES HAVE A RELATIONSHIP IN THEIR CORRECTIONS.   "
    "   NOTE: SINTHE AND COSTHE CAN BE CHANGED OUTSIDE THROUGH COMMON. "
    "   LVL IS A PARAMETER TELLING WHICH PARTICLES TO WORK WITH.       "
    "   THETA (SINTHE AND COSTHE) ARE ALWAYS RELATIVE TO THE DIRECTION "
    "   OF THE INCIDENT PARTICLE BEFORE ITS DIRECTION WAS ADJUSTED.    "
    "   THUS WHEN TWO PARTICLES NEED TO HAVE THEIR DIRECTIONS COMPUTED,"
    "   THE ORIGINAL INCIDENT DIRECTION IS SAVED IN THE VARIABLE A,B,C "
    "   SO THAT IT CAN BE USED ON BOTH CALLS."

    "   LVL=1 -- OLD PARTICLE, SAVE ITS DIRECTION AND ADJUST IT"
    "   LVL=2 -- NEW PARTICLE. ADJUST DIRECTION USING SAVED A,B,C"
    "   LVL=3 -- BREMSSTRAHLUNG GAMMA.  SAVE ELECTRON DIRECTION (NEXT  "
    "   TO TOP OF STACK), AND THEN ADJUST GAMMA DIRECTION."

    :NRK:
    GO TO (:OLD-PARTICLE:,:NEW-PARTICLE:,:BREMS-GAMMA:),LVL;
    "   LVL OUT-OF-BOUNDS IF HERE"   GO TO :ERROR:;

    :OLD-PARTICLE:
    A=U(NP);B=V(NP);C=W(NP);
    GO TO :ADJUST:;

    :BREMS-GAMMA:
    A=U(NP-1);B=V(NP-1);C=W(NP-1);

    :NEW-PARTICLE:
    $TRANSFER PROPERTIES TO (NP) FROM (NP-1);

    "   SEE H.H. NAGEL DISSERTATION FOR COORDINATE SYSTEM DESCRIPTION. "
    "   A ROTATION IS PERFORMED TO TRANSFORM DIRECTION COSINES OF THE  "
    "   PARTICLE BACK TO THE PHYSICAL FRAME (FROM THE TRANSPORT FRAME) "

    :ADJUST:
    SINPS2=A*A+B*B;
    "   If SINPS2 is small, no rotation is needed    "
    IF (SINPS2.LT.1.0E-20)["small polar angle case"
    U(NP)=SINTHE*COSPHI;
    V(NP)=SINTHE*SINPHI;
    W(NP)=C*COSTHE;    "fixed March 2001 from =COSTHE"
    ] "end small polar angle case"
    ELSE["large polar angle case"
    SINPSI=SQRT(SINPS2);
    US=SINTHE*COSPHI;
    VS=SINTHE*SINPHI;
    SINDEL=B/SINPSI;
    COSDEL=A/SINPSI;
    U(NP)=C*COSDEL*US-SINDEL*VS+A*COSTHE;
    V(NP)=C*SINDEL*US+COSDEL*VS+B*COSTHE;
    W(NP)=-SINPSI*US+C*COSTHE;
    ]"end large polar angle case"

    RETURN;


