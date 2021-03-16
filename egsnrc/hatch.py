import numpy
from egsnrc.params import *
from math import pi

# EMPTY CALLBACKS ----
hatch_user_input_init = None
need_photonuc_data = None
need_rayleigh_data = None

# ******************************************************************
#                                National Research Council of Canada
def hatch():
    """Load data for EGSnrc

    Setup which the user is expected to do before calling HATCH is:
    1. Set 'nmed' to the number of media to be used.
    2. Set the array 'media', which contains the names of the
        media that are desired.  the character format is a1, so
        that media(ib,im) contains the ib'th byte of the name of
        the im'th medium in a1 format.
    3. Set 'dunit', the distance unit to be used.
        dunit.gt.0 means value of dunit is length of distance unit
        centimeters.  dunit.lt.0 means use the radiation length of
        the abs(dunit)'th medium for the distance unit.
    4. Fill the array 'med' with the medium indices for the
        regions.
    5. Fill arrays 'ecut' and 'pcut' with the electron and photon
        cut-off energies for each region respectively.  setup will
        raise these if necessary to make them at least as large as
        the region's medium's ae and ap respectively.
    6. Fill 'med' array.  med(ir) is the medium index for region
        ir.  a zero medium index means the region is in a vacuum.
    7. Fill the array 'iraylr' with 1 for each region in which
        rayleigh (coherent) scattering is to be included.

    """

    # $ comin_hatch # DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
                    # COMIN/DEBUG,BOUNDS,BREMPR,EGS-VARIANCE-REDUCTION,
                        # ELECIN,MEDIA,MISC,PHOTIN,STACK,
                        # THRESH,UPHIIN,UPHIOT,USEFUL,USER,RANDOM/

    # Local HATCH variables in alphabetical order
    # $TYPE MBUF(72),MDLABL(8);
    # $REAL
    #     ACD   , "used to test goodness of sine-table look-up"
    #     ADEV  , "absolute deviation in sine-table look-up"
    #     ASD   , "used to test goodness of sine-table look-up"
    #     COST  , "cos(theta) from instrinsic library function"
    #     CTHET , "use to calculate cos(theta) according to look-up tables"
    #     DEL   , "leat squares delta for sine-table look-up"
    #     DFACT , "converts rl to dunits"
    #     DFACTI, "converts rl**-1 to dunits**-1"
    #     DUNITO, "units scaling varable"
    #     DUNITR, "saved value of dunit"
    #     FNSSS , "real form of integer nsinss"
    #     P     , "counter used in the pwr2i(i) = 1/2**(i - 1) construction"
    #     PZNORM, "used in $INITIALIZE-BREMS-ANGLE"
    #     RDEV  , "relative deviation in sine-table look-up"
    #     S2C2  , "sinthe**2 + costhe**2, used to test look-up table"
    #     S2C2MN, "min(s2c2)"
    #     S2C2MX, "max(s2c2)"
    #     SINT  , "sin(theta) from instrinsic library function"
    #     SX    , "sum of angles for least squared analysis of look-up table errors"
    #     SXX   , "sum**2 of angles for least square analysis of look-up table errors"
    #     SXY   , "sum of angle*sin(angle) for least squared analysis of look-up"
    #             "table errors"
    #     SY    , "sum of sin(angle) for least squared analysis of look-up table "
    #             "errors"
    #     WID   , "width of sine-table mesh sub-interval (sine-table algorithm)"
    #     XS    , "angle value in a sub-sub-interval (sine-table algorithm)"
    #     XS0   , "lower limit of a sub-sub-interval (sine-table algorithm)"
    #     XS1   , "upwer limit of a sub-sub-interval (sine-table algorithm)"
    #     XSI   , "beginning angle of a sun-interval (sine-table algorithm)"
    #     WSS   , "width of a sub-sub-interval (sine-table algorithm)"
    #     YS    , "sin(angle) for least squared analysis of look-up table errors"
    #     ZEROS(3); "zeros of sine, 0,pi,twopi"

    # $INTEGER
    #     I     , "generic do-loop variable"
    #     I1ST  , "flag = 0 on first pass"
    #     IB    , "do-loop variable used for reading the medium type"
    #     ID    , "integer value of -dunit, when dunit is negative"
    #     IE    , "do-loop variable for reading over elements in a compound/mixture"
    #     IL    , "do-loop variable used for reading the medium type"
    #     IM    , "do-loop variable looping over nmed, number of media"
    #     IRAYL , "Rayleigh switch read in from PEGS"
    #     IRN   , "do-loop variable over random set of sine-table look-ups"
    #     ISTEST, "flag that switches on test of sine function fit"
    #     ISUB  , "do-loop variable over sub-intervals of the sine look-up table"
    #     ISS   , "do-loop variable over sub-sub-intervals of the sine look-up table"
    #     IZ    , "used to locate an exact zero of a sun-interval mesh point in the"
    #             "sine-table look-up"
    #     IZZ   , "do-loop variable over the exact zeros of the sine-table look-up"
    #     J     , "do-loop variable looping over nmed, number of media"
    #     JR    , "do-loop variable looping over number of regions"
    #     LCTHET, "$SET INTERVAL index for cos(theta) from look-up table"
    #     LMDL  , "character width of medium header ' MEDIUM='"
    #     LMDN  , "character width of medium description"
    #     LTHETA, "$SET INTERVAL index for sin(theta) from look-up table"
    #     MD    , "temporary storage for the medium number"
    #     MXSINC, "number of intervals approximating the sine function"
    #     NCMFP , "array size input from PEGS. Currently 0, but probably intended"
    #             "to be cumulative electron mean free path. Presently unused."
    #     NEKE  , "array size input from PEGS."
    #             "Number of electron mapped energy intervals."
    #     NGE   , "array size input from PEGS."
    #             "Number of photon mapped energy intervals."
    #     NGRIM , "Rayleigh cross section array size."
    #     NISUB , "mxsinc - 2. Size of array with endpoints removed."
    #     NLEKE , "array size input from PEGS. Currently 0, but probably intended"
    #             "to be number of electron energy intervals below threshold."
    #             "Presently unused."
    #     NM    , "number of media found in the "
    #     NRANGE, "array size input from PEGS. Currently 0, but probably intended"
    #             "to be number of intervals in an array giving the electron range."
    #             "Presently unused."
    #     NRNA  , "number of random angles testing sine function fit"
    #     NSEKE , "array size input from PEGS. Currently 0, but probably intended"
    #             "to be number of electron small energy intervals. Presently unused."
    #     NSGE  , "array size input from PEGS. Currently 0, but probably intended"
    #             "to be number of gamma small energy intervals. Presently unused."
    #     NSINSS, "number of sub-intervals for each sine function interval"
    #     LOK($MXMED); "flag indicating that medium has been found in the PEGS "
    #                  "datafile"
    # character*256 tmp_string
    # ;integer*4  lnblnk1 #  in house lnblnk function becuase not all compilers
                        #  support this

    # DATA MDLABL/$S' MEDIUM='/,LMDL/8/,LMDN/24/,DUNITO/1./
    LMDL = 8
    LMDN = 24
    DUNITO = 1.0
    # DATA I1ST/1/,NSINSS/37/,MXSINC/MXSINC/,ISTEST/0/,NRNA/1000/
    I1ST = 1
    NSINSS = 37
    ISTEST = 0
    NRNA = 1000

    #    FORMAT STATEMENTS USED MULTIPLE TIMES IN SETUP
    # :INT:  FORMAT(1X,14I5)
    # :FLT:  FORMAT(1X,1PE14.5,4E14.5)
    # :BYTE: FORMAT(72A1)

    twopi = 2 * pi

    if i1st != 0:
        i1st = 0 # reset first time flag
        #    do first time initialization
        # --- Inline replace: $ HATCH_USER_INPUT_INIT; -----
        #  $ RNG-INITIALIZATION
        #  Have taken this out, (IK, Jan 2000). If the user does not initilize the
        #  rng before the first call to shower, the rng will initialize itself
        #  using the default seed and the default luxury level (which is defined
        #  via $ DEFAULT-LL).

        smaxir[smaxir <= 0] = 1e10
        # End inline replace: $ HATCH_USER_INPUT_INIT; ----

        # Now construct piecewise linear fit to sine function over the
        # interval (0, 5*pi/2).  Divide this interval into MXSINC sub-
        # intervals.  Each of these subintervals is then subdivided into
        # nsinss sub-sub-intervals.  The angles at the boundaries of
        # these sub-sub-intervals and their sines are used to compute
        # least squares coefficients for the subinterval.  An extra
        # subinterval on each side of the interval (0,5*pi/2) is included
        # for good measure.
        nisub = mxsinc-2
        fnsss = nsinss
        wid = pi5d2 / float(nisub)
        wss = wid/(fnsss-1.0)
        zeros = numpy.array(0, pi, twopi)

        for isub_m1 in range(MXSINC):  # Loop over subintervals
            isub = isub_m1 + 1
            # zero sums
            SX = 0.
            SY = 0.
            SXX = 0.
            SXY = 0.

            # Lower and upper limits
            xs0 = wid * float(isub-2)
            xs1 = xs0 + wid
            # Check to see if any zeros are in the interval
            iz = None
            for izz in range(3):
                if xs0 <= zeros[izz] <= xs1:
                    iz = izz
                    break

            if iz is None:
                xsi = xs0
            else:
                xsi = zeros[iz]

    #         DO ISS = 1,NSINSS [ # LOOP OVER SUB-SUBINTERVALS
    #             XS = WID*FLOAT(ISUB-2)+WSS*FLOAT(ISS-1)-XSI # ANGLE VALUE
    #             YS = SIN(XS+XSI) # SINE OF ANGLE
    #             SX = SX+XS # ACCUMULATE SUMS
    #             SY = SY+YS
    #             SXX = SXX+XS*XS
    #             SXY = SXY+XS*YS
    #         ] # END SUB-SUBINTERVAL LOOP

    #         # NOW COMPUTE LEAST SQUARES COEFFICIENTS
    #         if IZ != 0:
    #             [ # FORCE FIT THROUGH SINES' ZEROS,
    #             #           FOR SMALL REL.ERR.&GOOD
    #             # VALUES OF SINTHE/THETA NEAR ZERO
    #             SIN1(ISUB) = SXY/SXX
    #             SIN0(ISUB) = -SIN1(ISUB)*XSI
    #         else:
    #             DEL = FNSSS*SXX-SX*SX
    #             SIN1(ISUB) = (FNSSS*SXY-SY*SX)/DEL
    #             SIN0(ISUB) = (SY*SXX-SX*SXY)/DEL - SIN1(ISUB)*XSI

    #     ] # END SUB-INTERVAL LOOP

    #     SINC0 = 2.0  # SET COEFFICIENTS WHICH DETERMINE INTERVAL
    #     SINC1 = 1.0/WID

    #     # NOW TEST FIT, IF REQUESTED
    #     if ISTEST != 0:

    #         # FIRST TEST AT POINTS PREVIOUSLY COMPUTED, EXCLUDING
    #         # END SUBINTERVALS
    #         ADEV = 0.
    #         RDEV = 0.
    #         S2C2MN = 10.
    #         S2C2MX = 0.
    #         DO ISUB = 1,NISUB [
    #             DO ISS = 1,NSINSS [
    #             THETA = WID*FLOAT(ISUB-1)+WSS*FLOAT(ISS-1)
    #             CTHET = PI5D2-THETA


    #                 SINTHE =sin(THETA)
    #                 COSTHE =sin(CTHET)
    #             SINT = SIN(THETA)
    #             COST = COS(THETA)
    #             ASD = ABS(SINTHE-SINT)
    #             ACD = ABS(COSTHE-COST)
    #             ADEV = max(ADEV,ASD,ACD)
    #             if SINT != 0.0:

    #                 RDEV = max(RDEV,ASD/ABS(SINT))

    #             if COST != 0.0:

    #                 RDEV = max(RDEV,ACD/ABS(COST))

    #             S2C2 = SINTHE**2+COSTHE**2
    #             S2C2MN = min(S2C2MN,S2C2)
    #             S2C2MX = max(S2C2MX,S2C2)
    #             if ISUB < 11:

    #                 logger.info(4E20.7)',THETA,SINTHE,SINT,COSTHE,COST)


    #         ] # END OF FIXED INTERVAL TEST-OUTPUT RESULTS
    #         logger.info(2i5)',' SINE TESTS,MXSINC,NSINSS = ',MXSINC,NSINSS)
    #         $egs_info('(a,1PE16.8,3e16.8)',' ADEV,RDEV,S2C2(MN,MX) =',
    #                                         ADEV,RDEV,S2C2MN,S2C2MX)
    #         # NOW DO RANDOM TEST
    #         ADEV = 0.
    #         RDEV = 0.
    #         S2C2MN = 10.
    #         S2C2MX = 0.
    #         DO IRN = 1,NRNA [
    #             THETA = randomset()
    #             THETA = THETA*PI5D2
    #             CTHET = PI5D2-THETA
    #             SINTHE =sin(THETA)
    #             COSTHE =sin(CTHET)
    #             SINT = SIN(THETA)
    #             COST = COS(THETA)
    #             ASD = ABS(SINTHE-SINT)
    #             ACD = ABS(COSTHE-COST)
    #             ADEV = max(ADEV,ASD,ACD)
    #             if SINT != 0.0:

    #             RDEV = max(RDEV,ASD/ABS(SINT))

    #             if COST != 0.0:

    #             RDEV = max(RDEV,ACD/ABS(COST))

    #             S2C2 = SINTHE**2+COSTHE**2
    #             S2C2MN = min(S2C2MN,S2C2)
    #             S2C2MX = max(S2C2MX,S2C2)
    #         ] # END RANDOM ANGLE LOOP
    #         logger.info(i7,a)', ' TEST AT ',NRNA,' RANDOM ANGLES IN (0,5*PI/2)')
    #         $egs_info('(1PE16.8,3E16.8)',' ADEV,RDEV,S2C2(MN,MX) =',
    #                                     ADEV,RDEV,S2C2MN,S2C2MX)
    #     ] # END OF SINE TABLE TEST

    #     # NOW FILL IN POWER OF TWO TABLE.  PWR2I(I) = 1/2**(I-1)
    #     P = 1.
    #     DO I = 1,MXPWR2I [
    #         PWR2I(I) = P
    #         P = P/2.

    # ] # END OF FIRST TIME INITIALIZATION

    # # FILL IRAYLM ARRAY BASED ON IRAYLR INPUTS
    # # --- Inline replace: $ need_rayleigh_data; -----
    # if need_rayleigh_data:
    #     need_rayleigh_data()
    # else:

    #     DO J = 1,NMED [
    #     if goto_LOOP_OVER_REGIONS:  XXX  DO I = 1,MXREG [
    #     if IRAYLR(I) == 1.AND.MED(I).EQ.J:

    #     # REGION I = MEDIUM J AND WE WANT RAYLEIGH SCATTERING, SO
    #     # SET FLAG TO PICK UP DATA FOR MEDIUM J AND TRY NEXT MEDIUM.
    #     IRAYLM(J) = 1; EXIT :LOOP_OVER_REGIONS:;]
    #     # END OF REGION-LOOP]
    #     # END OF MEDIA-LOOP]
    # # End inline replace: $ need_rayleigh_data; ----

    # # Ali:photonuc, 2 lines
    # # FILL IPHOTONUCM ARRAY BASED ON IPHOTONUCR INPUTS
    # # --- Inline replace: $ need_photonuc_data; -----
    # if need_photonuc_data:
    #     need_photonuc_data()
    # else:

    #     IPHOTONUC = 0
    #     DO J = 1,NMED [
    #     if goto_LOOP_OVER_REGIONS_PHOTONUC:  XXX DO I = 1,MXREG [
    #     if IPHOTONUCR(I) == 1.AND.MED(I).EQ.J:

    #     # REGION I = MEDIUM J AND WE WANT PHOTONUCLEAR, SO
    #     # SET FLAG TO PICK UP DATA FOR MEDIUM J AND TRY NEXT MEDIUM.
    #     IPHOTONUCM(J) = 1; IPHOTONUC = 1; EXIT :LOOP_OVER_REGIONS_PHOTONUC:;]
    #     # END OF REGION-LOOP]
    #     # END OF MEDIA-LOOP]
    # # End inline replace: $ need_photonuc_data; ----
    # logger.info(i3)',' == = > Photonuclear flag: ', iphotonuc)

    # # NOW SEARCH FILE FOR DATA FOR REQUESTED MATERIALS
    # if ~is_pegsless:

    # REWIND KMPI

    # # explicit file name for HP compiler  Nov 23, 1996   DR
    # IUECHO = KMPO
    # NM = 0 # NUMBER OF MEDIA FOUND
    # DO IM = 1,NMED [
    # LOK(IM) = 0 # SET FLAG TELLING WHICH MEDIA ARE OK
    # # NOW TELL USER IF RAYLEIGH OPTION HAS BEEN REQUESTED
    # if IRAYLM(IM) == 1:
    #     logger.info(i3/)', ' RAYLEIGH OPTION REQUESTED FOR MEDIUM NUMBER',IM)


    # # Ali:photonuc, 1 block
    # DO IM = 1,NMED [
    # # TELL USER IF PHOTONUC HAS BEEN REQUESTED
    # if IPHOTONUCM(IM) == 1:

    #     logger.info(i3/)', ' PHOTONUCLEAR REQUESTED FOR MEDIUM NUMBER',IM)


    # if ~is_pegsless:

    # if goto_MEDIUM:  XXX
    # LOOP [# MEDIUM SEARCH LOOP

    # :MDLOOK:
    # LOOP [ # MEDIUM HEADER SEARCH LOOP
    #     # FIRST LOOK FOR MEDIUM HEADER
    #     READ(KMPI,:BYTE:,END = :MDNOMORE:)MBUF
    #     DO IB = 1,LMDL [
    #         if MBUF(IB) != MDLABL(IB):
    #             NEXT:MDLOOK:


    #     # HEADER MATCHES. NOW SEE IF IT IS ONE OF REQUESTED MEDIA
    #     :MDNAME:
    #     DO IM = 1,NMED [
    #         DO IB = 1,LMDN [
    #             IL = LMDL+IB
    #             if MBUF(IL) != MEDIA(IB,IM):

    #             NEXT:MDNAME:

    #             if IB == LMDN:

    #             EXIT:MDLOOK:


    #     ] # END :MDNAME: DO
    #     # NOT IN NAME TABLE, SO IGNORE IT
    # ] REPEAT # MDLOOK

    # # 'IM' IS THE INDEX OF THE MEDIUM READY TO BE READ
    # if LOK(IM) != 0:

    #     goto_MDLOOK = True
    #     break # XXX # WE ALREADY HAVE THIS ONE

    # LOK(IM) = 1
    # NM = NM+1 # SET FOUND FLAG AND STEP MEDIUM COUNTER

    # # NOW READY TO READ IN DATA FOR THIS MEDIUM
    # # $ UOUTPUT(KMPO)IM,MBUF;(' DATA FOR MEDIUM #',I3,', WHICH IS:',72A1)

    # # NOW PUT OUT LINES SHOWING COMPOSITION OF MEDIUM
    # # THE FOLLOWING LINE WAS CHANGED TO STORE THE ELEMENTAL COMPOSITION AFB 88/05/31
    # # $ UINPUT(KMPI)(MBUF(I),I = 1,5),RHO(IM),NE
    # # The next two lines were line prior to Dec 89 mods to get IUNRST
    # # $ UINPUT(KMPI)(MBUF(I),I=1,5),RHO(IM),NNE(IM)
    # # (5A1,5X,F11.0,4X,I2)
    # # following used to pick up IUNRST, IAPRIM and EPSTFL
    # # Problem is that GASP may or may not be printed, so we make
    # # a kludge which will work with all old data files
    # # FIRST WE ASSUME THERE IS NO GASP VALUE IN THE LINE
    # # Note that this reading scheme counts on there being an
    # # error when GASP does exist on the line--an error does
    # # occur on most compilers, however, we have found that on
    # # the rs6000 an error does not occur.  Instead, a warning
    # # is printed out and IUNRST,EPSTFL and IAPRIM are set to 0.
    # # This will make no difference in simulations but will cause
    # # a problem when running EXAMIN

    # #  IK: backspace(kmpi) fails under windows using g77 with I/O error
    # #   therefore we read the line in a temporary string and then
    # #   use memoty I/O to try to read with and without gasp there.

    # read(kmpi,'(a)',err=:hatch_read_error1:) tmp_string
    # goto_no_hatch_read_error1 = True
    # break # XXX
    # :hatch_read_error1:
    # logging.critical('***************** Error: ')

    #         logging.critical(''Error while reading pegs4 file'')

    #         logging.critical('***************** Quitting now.')

    #         sys.exit(1)


    # :no_hatch_read_error1:
    # read(tmp_string,1,ERR=:GASP_THERE:)
    # # READ(KMPI,1,ERR=:GASP_THERE:)
    # (MBUF(I),I=1,5),RHO(IM),NNE(IM),IUNRST(IM),EPSTFL(IM),IAPRIM(IM)
    # 1   FORMAT(5A1,5X,F11.0,4X,I2,9X,I1,9X,I1,9X,I1)
    # # IUNRST, EPSTFL AND IAPRIM ARE STORED IN COMIN ELECIN
    # goto_GASP_NOT_THERE = True
    # break # XXX

    # :GASP_THERE:
    # # WE MUST REREAD THE LINE WITH THE CORRECT FORMAT
    # # BACKSPACE(KMPI) # THIS BACKS UP ONE RECORD TO RE-READ IT
    # # READ(KMPI,2)

    # # The following output is only there because without it
    # # code compiled with the new gfortran GNU compiler
    # # fails with run time error. Another bug in their
    # # pre-alpha quality I/O system ----IK, Oct 26 2005
    # # write(6,*) 'Found medium with gas pressure'
    # logger.info('Found medium with gas pressure')
    # read(tmp_string,2)
    # (MBUF(I),I=1,5),RHO(IM),NNE(IM),IUNRST(IM),EPSTFL(IM),
    # IAPRIM(IM)
    # 2     FORMAT(5A1,5X,F11.0,4X,I2,26X,I1,9X,I1,9X,I1)

    # :GASP_NOT_THERE:

    # # THE FOLLOWING LINE WAS CHANGED AS WELL AFB 88/05/31
    # # $ UOUTPUT(KMPO)(MBUF(I),I=1,5),RHO(IM),NE
    # # ;$ UOUTPUT(KMPO)(MBUF(I),I=1,5),RHO(IM),NNE(IM)
    # # (5A1,',RHO=',1PG11.4,',NE=',I2,',COMPOSITION IS :')
    # # THE FOLLOWING LINE WAS CHANGED AS WELL AFB 88/05/31
    # # DO IE=1,NE[
    # DO IE = 1,NNE(IM)[
    #     # THE FOLLOWING LINE, COMMENTED OUT, WAS THE OLD WAY OF READING IN
    #     # THE ELEMENTAL COMPOSITION OF EACH MEDIUM. THE INFORMATION WAS NOT
    #     # PASSED ON TO EGS. IN THE PRESENT VERSION IT IS READ IN AND STORED
    #     # IN COMMON BREMPR. AFB 88/05/31.
    #     # READ(KMPI,:BYTE:)MBUF;WRITE(KMPO,:BYTE:)MBUF
    #     $UINPUT(KMPI)
    #     (MBUF(I),I=1,6),(ASYM(IM,IE,I),I=1,2),
    #     ZELEM(IM,IE),WA(IM,IE),PZ(IM,IE),RHOZ(IM,IE)
    #     (6A1,2A1,3X,F3.0,3X,F9.0,4X,F12.0,6X,F12.0)
    #     # $ UOUTPUT(KMPO)
    #     # (MBUF(I),I=1,6),(ASYM(IM,IE,I),I=1,2),
    #     # ZELEM(IM,IE),WA(IM,IE),PZ(IM,IE),RHOZ(IM,IE)
    #     # (6A1,2A1,',Z=',F3.0,',A=',F9.3,',PZ=',1PE12.5,',RHOZ=',1PE12.5)

    # # MEDIA AND THRESH
    # # $ ECHO  READ(KMPI,:FLT:) $ LGN(RLC,AE,AP,UE,UP(IM))
    # TE(IM) = AE(IM)-RM
    # THMOLL(IM) = TE(IM)*2. + RM

    # # ACTUAL ARRAY SIZES FROM PEGS
    # # $ ECHO  READ(KMPI,:INT:)
    # MSGE,MGE,MSEKE,MEKE,MLEKE,MCMFP,MRANGE(IM),IRAYL
    # NSGE = MSGE(IM)
    # NGE = MGE(IM)
    # NSEKE = MSEKE(IM)
    # NEKE = MEKE(IM)
    # NLEKE = MLEKE(IM)
    # NCMFP = MCMFP(IM)
    # NRANGE = MRANGE(IM)

    # # BREMPR
    # # $ ECHO  READ(KMPI,:FLT:)($ LGN(DL(I,IM)/1,2,3,4,5,6/),I=1,6)
    # # $ ECHO  READ(KMPI,:FLT:)DELCM(IM),($ LGN(ALPHI,BPAR,
    #     DELPOS(I,IM)),I=1,2)

    # # ELECIN
    # # $ ECHO  READ(KMPI,:FLT:)$ LGN(XR0,TEFF0,BLCC,XCC(IM))
    # # $ ECHO  READ(KMPI,:FLT:)$ LGN(EKE(IM)/0,1/)
    # # $ ECHO  READ(KMPI,:FLT:)
    # ($LGN(ESIG,PSIG,EDEDX,PDEDX,EBR1,PBR1,PBR2,
    #     TMXS(I,IM)/0,1/),I=1,NEKE)

    # # PHOTIN
    # # $ ECHO  READ(KMPI,:FLT:)EBINDA(IM),$ LGN(GE(IM)/0,1/)
    # # $ ECHO  READ(KMPI,:FLT:)($ LGN(GMFP,GBR1,GBR2(I,IM)/0,1/),I=1,NGE)

    # # PHOTIN (CONTINUED)---OPTIONAL RAYLEIGH SCATTERING INPUT

    # /* Leave this for compatibility with existing pegs4 data sets.  */
    # if IRAYL == 1:

    #     # $ ECHO  READ(KMPI,:INT:) NGR(IM)
    #     NGRIM = NGR(IM)
    #     # $ ECHO  READ(KMPI,:FLT:)$ LGN(RCO(IM)/0,1/)
    #     # $ ECHO  READ(KMPI,:FLT:)($ LGN(RSCT(I,IM)/0,1/),I=1,NGRIM)
    #     # $ ECHO  READ(KMPI,:FLT:)($ LGN(COHE(I,IM)/0,1/),I=1,NGE)
    #     # if(IRAYLM(IM) != 1) [
    #     $egs_info('(a,i3,a)', ' Rayleigh data available for medium',
    #             IM, ' in PEGS4 data set.')
    #     # ]

    # /*******************************************************************
    # Rayleigh data picked up directly from pgs4form.data or user-supplied
    # ff file in egs_init_rayleigh unless user wants to use PEGS4 data.
    # *********************************************************************/
    # if IRAYLM(IM) == 1:
    #         [ # Rayleigh data requested for medium IM
    #     if IRAYL != 1:
    #         [ # No data in PEGS4
    #         if toUpper(photon_xsections(:lnblnk1(photon_xsections))) == 'PEGS4':
    #             [# Rayleigh not possible
    #             $egs_fatal('(a,i3 /,a /,a)',
    #             ' IN HATCH: REQUESTED RAYLEIGH OPTION FOR MEDIUM',
    #             IM,' BUT RAYLEIGH DATA NOT INCLUDED IN PEGS4 FILE.',
    #             ' YOU WILL NOT BE ABLE TO USE THE PEGS4 DATA WITH RAYLEIGH ON!')
    #         else:
    #         $egs_warning('(a,i3 /,a)',
    #         ' IN HATCH: REQUESTED RAYLEIGH OPTION FOR MEDIUM',
    #         IM,' BUT RAYLEIGH DATA NOT INCLUDED IN PEGS4 FILE.')

    #     else:
    #         if toUpper(photon_xsections(:lnblnk1(photon_xsections))) == 'PEGS4':
    #             [ # PEGS4 data selected
    #             # ***********************************************************
    #             # Preparing data for new Rayleigh angular sampling when using
    #             # the pegs4 data set,
    #             # ***********************************************************
    #             call egs_init_rayleigh_sampling(IM)

    #         # ELSE[Taking photon data from either si,epdl,xcom or user]


    # /*******************************************************************/

    # # THAT'S ALL FOR THIS MEDIUM
    # ] UNTIL NM.GE.NMED # LOOP UNTIL WE HAVE ENOUGH.  END :MEDIUM: LOOP

    # CLOSE (UNIT=KMPI)

    # # WE NOW HAVE DATA FOR ALL MEDIA REQUESTED.  NOW DO DISTANCE UNIT
    # # CHANGE.  DATA FROM PEGS IS IN UNITS OF RADIATION LENGTHS.
    # # EGS IS RUN IN UNITS OF 'DUNIT' CENTIMETERS, if DUNIT > 0
    # # OR IN UNITS OF RLC(-DUNIT) CENTIMETERS if DUNIT < 0.
    # # THAT IS, A NEGATIVE DUNIT MEANS UNIT IS TO BE THE RADIATION
    # # LENGTH OF THE MEDIUM WHOSE INDEX IS -DUNIT
    # DUNITR = DUNIT # SAVE REQUESTED
    # if DUNIT < 0.0:

    # ID = MAX0(1,MIN0(MXMED,int(-DUNIT)))
    # DUNIT = RLC(ID)

    # if DUNIT != 1.0:

    # $egs_info('(a,1PE14.5,E14.5,a)',' DUNIT REQUESTED&USED ARE: ',
    #         DUNITR,DUNIT,'(CM.)' )

    # DO IM = 1,NMED [
    # DFACT = RLC(IM)/DUNIT # CONVERTS RL TO DUNITS
    # DFACTI = 1.0/DFACT # CONVERT RL**-1 TO DUNITS**-1

    # FOR I = 1 TO MEKE(IM) [
    # ESIG= $LGN0,1(ESIG*DFACTI;PSIG = PSIG*DFACTI;EDEDX = EDEDX*DFACTI;PDEDX(I = PDEDX(I*DFACTI;IM) = IM)*DFACTI
    # TMXS0,1(I= $LGN(TMXS(I*DFACT;IM) = IM)*DFACT

    # TEFF0(IM)= TEFF0(IM)*DFACT
    # BLCC(IM)= BLCC(IM)*DFACTI
    # XCC(IM)= XCC(IM)*SQRT(DFACTI)
    # RLDU(IM) = RLC(IM)/DUNIT
    # FOR I = 1 TO MGE(IM) [
    # GMFP0,1(I= $LGN(GMFP(I*DFACT;IM)=IM)*DFACT

    # ] # END IM DO

    # # SCALE VACDST.  UNDO PREVIOUS SCALE, THEN DO NEW.
    # VACDST = VACDST*DUNITO/DUNIT
    # DUNITO = DUNIT # SAVE OLD DUNIT

    # ] # end regular pegs4 intake
    # else:

    # logger.info(' PEGSLESS INPUT.  CALCULATING ELECTRON CROSS-SECTIONS.')


    # $egs_fatal('(a/a)',' Code cannot be run in pegsless mode.',
    # ' Compile with required files and try again.')

    # # NOW MAKE SURE ECUT AND PCUT ARE NOT LOWER THAN ANY AE OR AP
    # # ALSO SET DEFAULT DENSITIES

    # DO JR = 1,MXREG [
    #     MD = MED(JR)
    #     if (MD >= 1).AND.(MD <= NMED):
    #         [# IT IS LEGAL NON-VACUUM MEDIUM.
    #         ECUT(JR) = max(ECUT(JR),AE(MD))
    #         PCUT(JR) = max(PCUT(JR),AP(MD))
    #         # USE STANDARD DENSITY FOR REGIONS NOT SPECIALLY SET UP
    #         if RHOR(JR) == 0.0)[RHOR(JR) = RHO(MD:
    #             ;]


    # # BREMSSTRAHLUNG ANGULAR DISTRIBUTION INITIALIZATION - DEFAULT IS NULL
    # # NEXT LINE ADDED AFB 88/05/31

    # ; if(IBRDST == 1)[
    #         DO IM = 1,NMED[
    #             ZBRANG(IM) = 0.0;
    #             PZNORM = 0.0
    #             DO IE = 1,NNE(IM)[
    #                 ZBRANG(IM) =
    #                 ZBRANG(IM)+PZ(IM,IE)*ZELEM(IM,IE)*(ZELEM(IM,IE)+1.0)
    #                 PZNORM = PZNORM+PZ(IM,IE)

    #                 ZBRANG(IM) = (8.116224E-05)*(ZBRANG(IM)/PZNORM)**(1./3.)
    #                 LZBRANG(IM) = -log(ZBRANG(IM))


    # # PAIR ANGULAR DISTRIBUTION INITIALIZATION - DEFAULT IS NULL
    # # NEXT LINE ADDED AFB 91/05/29

    # ;    if(IPRDST > 0)[
    #         DO IM = 1,NMED[
    #             ZBRANG(IM) = 0.0;PZNORM = 0.0
    #             DO IE = 1,NNE(IM)[
    #                 ZBRANG(IM) =
    #                 ZBRANG(IM)+PZ(IM,IE)*ZELEM(IM,IE)*(ZELEM(IM,IE)+1.0)
    #                 PZNORM = PZNORM+PZ(IM,IE)

    #                 ZBRANG(IM) = (8.116224E-05)*(ZBRANG(IM)/PZNORM)**(1./3.)


    # #  See if user has requested PEGS4 photon cross section data
    # if toUpper(photon_xsections(:lnblnk1(photon_xsections))) == 'PEGS4':

    # $egs_warning('(6(a/))','Using photon data from PEGS4 file!!!',
    # 'However, the new Rayleigh angular sampling will be used.',
    # 'The original EGS4 angular sampling undersamples large scattering ',
    # 'angles. This may have little impact as Rayleigh scattering ',
    # 'is forward peaked.',
    # '*********************************************************')

    # else:
    # # Ali:photonuc, 2 lines
    #     call egs_init_user_photon(photon_xsections,comp_xsections,
    #     photonuc_xsections,xsec_out)
    # #  call egs_init_user_photon(photon_xsections,comp_xsections,xsec_out)

    # call mscati # Initialize new MS, step-sizes, etc, IK Oct 97

    # # Calling order of the subroutines below is important when using
    # # detailed atomic relaxation in order to use the binding energies
    # # corresponding to the requested photon cross section library
    # if  eadl_relax and photon_xsections == 'xcom' :

    # call init_compton # Initialize bound Compton scattering
    # call EDGSET(1,1) # Initialize relaxations and photo-absorption data
    # else:
    # call EDGSET(1,1) # Initialize relaxations and photo-absorption data,
    #                     # if requested
    # call init_compton # Initialize bound compton scattering, IK, Jan 99
    #                     # if requested

    # if  xsec_out == 1 and eadl_relax:

    # call egs_print_binding_energies

    # call fix_brems # Re-calculate dl1,... for the different technique
    #                 # employed in BREMS. Note that the old EGS sampling
    #                 # technique for BREMS had a bug that shows up only
    #                 # if AP is not much smaller than electron kinetic energy

    # if  ibr_nist >= 1 :

    # call init_nist_brems

    #                 # initializes the sampling tables and modifies the total
    #                 # brems cross sections if the NIST brems data base is to
    #                 # be used

    # if  pair_nrc == 1 :

    # call init_nrc_pair

    # #  Load and initialize EII data if needed.
    # call eii_init

    # #  Load and initialize the triplet data if needed
    # call init_triplet

    # # SETUP IS NOW COMPLETE
    # if NMED == 1:

    # logger.info('EGSnrc SUCCESSFULLY ''HATCHED'' FOR ONE MEDIUM.')
    # else:
    # $egs_info('(a,i5,a)',
    #             'EGSnrc SUCCESSFULLY ''HATCHED'' FOR ',NMED,' MEDIA.')

    # RETURN

    # if goto_MDNOMORE:  XXX
    # $egs_info('(a,i2//,a/,a/)', ' END OF FILE ON UNIT ',KMPI,
    # ' PROGRAM STOPPED IN HATCH BECAUSE THE',
    # ' FOLLOWING NAMES WERE NOT RECOGNIZED:')
    # DO IM = 1,NMED [
    # if LOK(IM) != 1:

    #     logger.info(a,24a1,a)','''',(MEDIA(I,IM),I=1,LMDN),'''')


    # STOP
    # # END OF def HATCH   END:


def key_vals(line) -> dict:
    """Split a comma-separated line of key=value pairs into converted values

    Each value converted to `float` or `int`, if valid, else left as `str`
    """
    parts = [
        kv.strip() for kv in line.split(",")
        if '=' in kv
    ]
    di = {}
    for part in parts:
        k, v = part.split("=")
        if '.' in v:
            try:
                val = float(v)
            except ValueError:  # in case of string with '.'
                val = v
        else:
            try:
                val = int(v)
            except ValueError:
                try:
                    val = float(v)  # still try float in case of e.g. '1e8'
                except ValueError:
                    val = v  # string
        di[k.lower()] = val
    return di


def read_floats(f, num):
    """Read floats 5 per line until `num` read"""
    li = []
    while num > 0:
        floats = [float(x) for x in next(f).split()]
        expected = 5 if num >= 5 else num
        assert len(floats) == expected
        li.extend(floats)
        num -= expected

    return numpy.array(li)


def read_pegs(filename: str, requested_media: set):
    """Get the physics data for the requested mediums

    Return dict of medium name to its properties/data
    """

    # Start with all requested media, remove them as they are found
    needed_media = [x.lower() for x in requested_media]
    media_dict = {}
    with open(filename, 'r') as f:
        while needed_media:
            # Find ' MEDIUM=' line
            # e.g. '  MEDIUM=NAI                     ,STERNCID=NAI'
            line = next(f)
            while not line.startswith(" MEDIUM="):
                line = next(f)
            medium_dict = key_vals(line)
            # Check if don't need or already have (orig Mortran skipped if had)
            medium_name = medium_dict['medium'].lower()
            if medium_name not in needed_media:
                continue  # find next MEDIUM= in file

            # Found a medium we need, read all its data
            media_dict[medium_name] = medium_dict
            # First, med_type, rho, num_elems
            # e.g. ' COMP,RHO= 3.6700E+00,NE= 2'
            # First has no "=", add in to get 'type' in dict of parsed line
            medium_dict.update(key_vals("TYPE=" + next(f)))
            num_elems = medium_dict['ne']
            # Get each elem details - ASYM,Z,A,PZ, RHOZ
            # -> in Fortran arrays these are:
            # ASYM(IM,IE,I=1,2) (2chr symbl), ZELEM(IM,IE), WA(IM,IE),
            #     PZ(IM,IE), RHOZ(IM,IE), where IM = medium index, IE=elem index
            # e.g. for NaI:
            # ' ASYM=NA,Z=11.,A=   22.990,PZ= 1.00000E+00,RHOZ= 2.29898E+01
            # ' ASYM=I ,Z=53.,A=  126.904,PZ= 1.00000E+00,RHOZ= 1.26904E+02
            elem_details = [key_vals(next(f)) for i in range(num_elems)]
            medium_dict['elements'] = elem_details

            # Get medium data.
            # rlc = radiation length in cm
            # ae, ap = creation energy thresholds electron/photon
            # ue, up = upper energy in pegs dataset electron/photon
            # e.g. for NaI:
            # '    2.58633E+00   7.00000E-01   1.00000E-02   5.05110E+01   5.00000E+01'
            keys = ("rlc", "ae", "ap", "ue", "up")
            vals = [float(x) for x in next(f).split()]
            medium_dict.update(zip(keys, vals))

            # Get info on physics data
            # In Mortran $LGN(MSGE,MGE,MSEKE,MEKE,MLEKE,MCMFP,MRANGE(IM)),IRAYL
            # Variable explanations from Mortran COMMON/MEDIA in egsnrc.macros:
            # MSGE ?  MSEKE ?   MLEKE ?,   MCMFP ?,   MRANGE ?
            # MGE: num photon mapped energy intervals (?'Mapped Gamma Energies'?)
            # MEKE: num e mapped energy intervals (?'Mapped Electron KEs'?)
            # IRAYL:  Rayleigh switch from PEGS
            # e.g. for NaI in tutor_data:
            # '     0  199    0   46    0    0    0    1    0'
            # Note, is one extra int in vals, zip goes to shorter list
            keys = "msge mge mseke meke mleke mcmfp mrange irayl".split()
            vals = [int(x) for x in next(f).split()]
            medium_dict.update(zip(keys, vals))

            # BREMPR
            # Mortran: ($LGN(DL(I,IM)/1,2,3,4,5,6/),I=1,6)
            #     and DELCM(IM),($LGN(ALPHI,BPAR,DELPOS(I,IM)),I=1,2)
            floats = read_floats(f, 6*6)
            floats.shape = (6, 6)
            medium_dict.update({
                f"dl{i+1}": floats[:, i].transpose()
                for i in range(6)
            })
            floats = read_floats(f, 1 + 3*2)
            medium_dict.update({
                'delcm': floats[0],
                'alphi': numpy.array((floats[1], floats[4])),
                'bpar': numpy.array((floats[2], floats[5])),
                'delpos': numpy.array((floats[3], floats[6]))
            })

            # Now have all info on this medium, remove from list to find
            needed_media.remove(medium_name)

            # ELECIN
            # Mortran: $LGN(XR0,TEFF0,BLCC,XCC(IM))
            #     and  $LGN(EKE(IM)/0,1/)
            #     and  $LGN(ESIG,PSIG,EDEDX,PDEDX,EBR1,PBR1,PBR2,TMXS(I,IM)/0,1/),I=1,NEKE)
            keys = "xr0,teff0,blcc,xcc".split(",")
            floats = read_floats(f, 4)
            medium_dict.update(zip(keys,floats))

            keys = "eke0 eke1".split()
            floats = read_floats(f, 2)
            medium_dict.update(zip(keys,floats))

            neke = medium_dict['meke']
            floats = read_floats(f, 8*2*neke)
            floats.shape = (neke, 8, 2)

            arr0s = floats[:,:,0].transpose()
            arr1s = floats[:,:,1].transpose()
            var_names = "esig psig ededx pdedx ebr1 pbr1 pbr2 tmxs".split()
            medium_dict.update(
                {var_name+"0": arr0s[i] for i, var_name in enumerate(var_names)}
            )
            medium_dict.update(
                {var_name+"1": arr1s[i] for i, var_name in enumerate(var_names)}
            )

            # PHOTIN
            # EBINDA(IM),$LGN(GE(IM)/0,1/);
            # $LGN(GMFP,GBR1,GBR2(I,IM)/0,1/),I=1,NGE);
            nge = medium_dict['mge']
            floats = read_floats(f, 1+2)
            medium_dict.update({
                'ebinda': floats[0],
                'ge0': floats[1],
                'ge1': floats[2]
            })

            floats = read_floats(f, 3*2*nge)
            floats.shape = (nge, 3, 2)

            arr0s = floats[:,:,0].transpose()
            arr1s = floats[:,:,1].transpose()
            var_names = "gmfp gbr1 gbr2".split()
            medium_dict.update(
                {var_name+"0": arr0s[i] for i, var_name in enumerate(var_names)}
            )
            medium_dict.update(
                {var_name+"1": arr1s[i] for i, var_name in enumerate(var_names)}
            )

            if medium_dict['irayl'] != 1:
                continue

            # Rayleigh
            # READ(KMPI,:INT:) NGR(IM);
            ngr = int(next(f).strip())

            # READ(KMPI,:FLT:)$LGN(RCO(IM)/0,1/);
            medium_dict['rco0'], medium_dict['rco1'] = read_floats(f, 2)

            # READ(KMPI,:FLT:)($LGN(RSCT(I,IM)/0,1/),I=1,NGRIM);
            floats = read_floats(f, 2*ngr)
            floats.shape = (ngr, 2)
            medium_dict.update({
                'rsct0': floats[:, 0].transpose(),
                'rsct1': floats[:, 1].transpose(),
            })

            # READ(KMPI,:FLT:)($LGN(COHE(I,IM)/0,1/),I=1,NGE);
            floats = read_floats(f, 2*nge)
            floats.shape = (nge, 2)
            medium_dict.update({
                'cohe0': floats[:, 0].transpose(),
                'cohe1': floats[:, 1].transpose(),
            })

    return media_dict
