import logging
logger = logging.getLogger("egsnrc")

from egsnrc.commons import *


log_it = logger.debug  # can set this if wish to change to different log level

icount = 0
jhstry = 1
graph_unit = -1

_2X = "  "
_3X = "   "
_6X = "      "
_7X = "       "
header = "".join(
    ("\n\n"," "*39,' NP',_3X,'ENERGY  Q REGION    X',_7X,
    'Y',_7X,'Z',_6X,'U',_6X,'V',_6X,'W',_6X,'LATCH',_2X,'WEIGHT\n'
    )
)

# Simple iargs that just print a message and current state
arg_messages = {
    1: ' Discard  AE,AP<E<ECUT',
    2: ' Discard  E<AE,AP',
    3: ' Discard -user request',
    6: ' bremsstrahlung  about to occur',
    8: ' Moller   about to occur',
    10: ' Bhabba   about to occur',
    12: ' Positron about to decay in flight',
    15: ' Pair production about to occur',
    17: ' Compton  about to occur',
    19: ' Photoelectric about to occur',
    24: ' Rayleigh scattering occured',
    25: 'Fluorescent X-ray created',
    26: 'Coster-Kronig e- created',
    27: 'Auger electron created',
    28: ' Positron will annihilate at rest',
}


interaction_rejected_iargs = (9, 11, 13, 14, 16, 18)
rejected_message = '           Interaction rejected'
second_message = {
    9: '          Resulting electrons',
    11:'          Resulting e- or e+',
    13:'          Resulting photons',
    14: ' Positron annihilates at rest',
}


def std_data(msg, ip, ke):
    """ip must be already -1 for 0-based arrays in Python"""
    return (
        f"{msg:<35}:{ip:5}{ke:9.3f}{iq[ip]:4}{ir[ip]:4}"
        f"{x[ip]:8.3f}{y[ip]:8.3f}{z[ip]:8.3f}"
        f"{u[ip]:7.3f}{v[ip]:7.3f}{w[ip]:7.3f}"
        f"{latch[ip]:10}{wt[ip]:10.3E}"
    )


def watch(iarg, iwatch):
    """====================================================================

      A general purpose auxiliary routine for use with the EGSnrc system

      It prints out information about the particle transport

        For iwatch = 1 it prints information about each discrete interaction
        For iwatch = 2 or 3 it prints information about each step as well
        For iwatch = 4 it prints graphing data for use with EGS_Windows


     Routine is used via two mandatory and 1 optional call from the user's
           code

    1)The routine must be initialized by a call with iarg=-99 before the first
           call to SHOWER. It should be after all inputs are in place.
    2)The routine must be called near the beginning of the AUSGAB subroutine
           IF (iwatch > 0 ) CALL WATCH(iarg,iwatch)
    3)The routine may be called at the end of each history with iarg = - 1 so
           a message will get printed stated history is complete

     Since WATCH cannot output values related to the initial values in a
     shower call, it is useful to also put something like the following
     immediately prior to the CALL SHOWER stmt
            logger.debug(
            "\n INITIAL SHOWER VALUES             :"
            f"    1{ei:9.3f}{iqin:4}{irin:4}"
            f"{xin:8.3f}{yin:8.3f}{zin:8.3f}"
            f"{uin:7.3f}{vin:7.3f}{win:7.3f}"  # should be 8.3 like x,y,z but get extra spaces
            f"{latchi:10}{wtin:10.3E}",
            )

     Note ei is the kinetic energy of the incident particle


    The routine uses up to 132 columns for output.

      JAN 1984  GENERALIZED VERSION WITH INITIALIZATION
                               DAVE ROGERS NRCC
      JUN 1987  PUT IN iwatch = 4 OPTION     AFB
      JUL 1988  COMPATIBLE WITH X-RAY FLUORESCENCE  DWOR
      SEP 1990  ADDED ENERGY OUTPUT TO iwatch = 4 OPTION     AFB
      OCT 1990  UNIX compatible carriage control   DWOR
      JAN 2000  Rewritten to output relaxation particles and also
                so some of the output makes more sense BW
      FEB 2021  Transpiled to Python  DLM
    """

    global icount, jhstry, graph_unit

    ku = 13; kr = 0; ka = 1  # graph file params
    if iarg == -99:
        # Initialize flags so we will get calls thru AUSGAB
        iausfl[:] = 1
        iausfl[21:24] = 0

    if iarg == -1:
        # main is assumed to call AUSGAB with iarg=-1 at end of history
        if iwatch == 4:
            raise NotImplementedError("Python watch() does not handle graph output yet")
            # if  graph_unit < 0:
            #     graph_unit = egs_open_file(ku,kr,ka,'.egsgph')
            # WRITE(graph_unit,:GRAPHICS_FORMAT:) 0,0,0,0.0,0.0,0.0,0.0,jhstry
            # jhstry += 1
        else:
            log_it(f" END OF HISTORY {jhstry:8}   " + "*"*40)
            jhstry += 1
            icount += 2
            return

    if iwatch != 4 and (icount >= 50 or icount == 0 or iarg == -99) :
        # PRINT HEADER
        icount = 1
        log_it(header)

    if iwatch == 4 and iarg >= 0 and (iarg != 5):
        # GRAPHICS OUTPUT
        raise NotImplementedError("Python watch() does not handle graph output yet")
        # if  graph_unit < 0 ) graph_unit == egs_open_file(ku,kr,ka,'.egsgph':
        # WRITE(graph_unit,:GRAPHICS_FORMAT:) NP,iq[np_m1],ir[np_m1],x[np_m1],y[np_m1],z[np_m1],e[np_m1]
        # :GRAPHICS_FORMAT:FORMAT(2I4,1X,I6,4G15.8,I12)

    if iarg == 5 or iarg < 0:
        return

    if iwatch == 4:
        return # none of the rest needed for graphics output

    np_m1 = np - 1  # ** 0-based arrays
    npold_m1 = npold - 1
    ke = e[np_m1]
    if iq[np_m1] != 0:
        ke = e[np_m1] - prm

    if iarg == 0:
        if iwatch == 2:
            icount += 1
            log_it(std_data('           STEP ABOUT TO OCCUR', np_m1, ke))
        else:
            return

    if iarg in arg_messages:
        icount += 1
        log_it(std_data(arg_messages[iarg], np_m1, ke))
    elif iarg == 4:
        msg = '         Local energy deposition'
        # icount += 1  # not done in original EGSnrc, but should be?
        log_it(f"{msg:<35}:{edep:12.5f} MeV in region {ir[np_m1]:6}")
    elif iarg == 7:
        if nbr_split == 1:  # no splitting or SBS is on in BEAMnrc
            for ip_m1 in range(npold_m1, np):
                if iq[ip_m1] == -1:
                    ke = e[ip_m1] - rm
                    msg = '          Resulting electron'
                else:
                    ke = e[ip_m1]
                    msg = '          Resulting photon'
                icount += 1
                log_it(std_data(msg, ip_m1, ke))
        else:  # splitting case--e- is always at NPold
            ke = e[npold_m1] - rm
            icount += 1
            log_it(std_data('          Resulting electron', npold_m1, ke))
            for ip_m1 in range(npold+1-1, np):  # 0-based
                ke= e[ip_m1]
                # print label info for first one only"
                msg = '          Split photons' if ip_m1 == npold else ""
                icount += 1
                log_it(std_data(msg, ip_m1, ke))
    elif iarg in interaction_rejected_iargs:
        msg = rejected_message
        if np == npold and (i_survived_rr == 0 if iarg in (16,18) else True):
            icount += 1
            log_it(std_data(msg, np_m1, ke))
        elif iarg in (9, 11, 13, 14):
            for ip_m1 in range(npold_m1, np):
                ke = e[ip_m1] - abs(iq[np_m1]) * rm
                msg = second_message[iarg] if ip_m1 == npold_m1 else ""
                icount += 1
                log_it(std_data(msg, ip_m1, ke))
        elif iarg == 16:  # pair production
            if np == npold and i_survived_rr > 0:  # we have cleared the stack
                log_it(
                    f'          Russian Roulette eliminated {i_survived_rr:2}'
                    f' particle(s) with probability {prob_rr:8.5f}'
                )
                icount += 1
                log_it(std_data('Now on top of stack', np_m1, ke))
            else:
                for ip_m1 in range(npold-1,np):
                    ke = e[ip_m1] - abs(iq[ip_m1]) * rm
                    msg = 'resulting pair' if ip_m1 == npold-1 else ""
                    icount += 1
                    log_it(std_data(msg, ip_m1, ke))
        elif iarg == 18:
            if np > npold:  # have not cleared the stack with rus rou
                for ip_m1 in range(npold-1,npold+1):
                    ke = e[ip_m1] - abs(iq[ip_m1]) * rm
                    msg = 'compton electron created' if iq[ip_m1] != 0 else 'compton scattered photon'
                    icount += 1
                    log_it(std_data(msg, ip_m1, ke))

            if(i_survived_rr > 0):  # whether the stack has been cleared or not
                log_it(
                    f'          Russian Roulette eliminated {i_survived_rr:2}'
                    f' particle(s) with probability {prob_rr:8.5f}'
                )
                icount += 1
                log_it(std_data('Now on top of stack', np_m1, ke))

        else:
            # XXX defensive in case of log errors, and `else` above should work
            raise ValueError("Should only have been iarg 16 or 18")
    elif iarg == 20:
        if npold == np and iq[np_m1] == 0 and i_survived_rr == 0:
            msg = (
                '          Photon energy below N-shell\n'
                '          Photon discarded'
            )
            icount += 1
            log_it(std_data(msg, np_m1, ke))
        elif iq[npold_m1] == -1 and i_survived_rr == 0:
            ke= e[npold_m1] - rm
            msg = '         Resulting photoelectron'
            icount += 1
            log_it(std_data(msg, npold_m1, ke))
        elif i_survived_rr > 0:  # done some russian roulette
            if np == npold-1 or iq[npold_m1] != -1:
                if i_survived_rr > 1:  # eliminated more than the photoelectron
                    log_it(
                        f'          Russian Roulette eliminated {i_survived_rr - 1:2}'
                        f' particle(s) with probability {prob_rr:8.5f}'
                    )
                log_it(
                    f'          Russian Roulette eliminated resulting photoelectron'
                    f' with probability {prob_rr:8.5f}'
                )
            else:  # NPold could hold the photoelectron
                ke = e[npold_m1] - rm
                msg = '         Resulting photoelectron?'
                icount += 1
                log_it(std_data(msg, npold_m1, ke))
                log_it(
                    f'          Russian Roulette eliminated {i_survived_rr:2}'
                    f' particle(s) with probability {prob_rr:8.5f}'
                )
            msg = '         Now on top of stack'
            icount += 1
            log_it(std_data(msg, np_m1, ke))

    if iarg == 0 and iwatch == 2:
        msg = '    USTEP,TUSTEP,VSTEP,TVSTEP,EDEP'
        log_it(
            f"{msg:<35}:    {ustep:13.4E},{tustep:13.4E},{vstep:13.4E},"
            f"{tvstep:13.4E},{edep:13.4E}"
        )
        icount += 1

    if np == 1 or iarg == 0:
        return

    if iarg <= 3:
        N = np-1
        N_m1 = N - 1
        ke = e[N_m1] - abs(iq[N_m1]) * rm
        msg = '         Now on top of stack'
        icount += 1
        log_it(std_data(msg, N_m1, ke))
