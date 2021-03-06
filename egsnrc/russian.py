from .commons import *
from .params import *


# This macro implements Russian Roulette (most useful  with brems splitting)
# It is more efficient than having the user do it via AUSGAB since it avoids
# considerable handling of the particles by ELECTR
# The user must set i_play_RR (defaults to 0) and prob_rr
# Both are in COMIN EGS-VARIANCE-REDUCTION
#
# Note that this macro is called as $ PLAY RUSSIAN ROULETTE WITH ELECTRONS...
# Note also that subroutine pair has its own, internal version

def russian_roulette(ip):
    """Implements Russian Roulette for electrons

    (most useful  with brems splitting)
    This is more efficient than having the user do it via AUSGAB since it avoids
    considerable handling of the particles by ELECTR
    The user must set `i_play_RR` (defaults to 0) and `prob_rr`
    Both are in egs_variance_reduction

    Note also that subroutine pair has its own, internal version

    Parameters
    ----------
    ip: int
        Index of particle
    """
    ip_m1 = ip - 1  # ** 0-based
    i_survived_rr = 0 # flag all survive
    if i_play_rr == 1:
        if prob_rr <= 0:
            if n_rr_warning < MAX_RR_WARNING:
                n_rr_warning = n_rr_warning + 1
                logger.warning(
                    '**** Warning, attempt to play Roussian Roulette with'
                    f' prob_rr<=0! {prob_rr:14.6f}'
                )
            return

        while True:  # handle all particles from p1 to np
            if iq[ip_m1] != 0:
                # i.e. charged particles
                if randomset() < prob_rr:
                    # particle survives
                    wt[ip_m1] /= prob_rr
                    ip += 1 # increase local pointer
                    ip_m1 += 1
                else:
                    egs_vr.i_survived_rr += 1
                    if ip < np:
                        # => replace it with last particle on stack
                        e[ip_m1] = e[np_m1]
                        iq[ip_m1] = iq[np_m1]
                        wt[ip_m1] = wt[np_m1]
                        u[ip_m1] = u[np_m1]
                        v[ip_m1] = v[np_m1]
                        w[ip_m1] = w[np_m1]

                    stack.np -= 1 # reduce stack by one=> particle gone
                    np_m1 -= 1
                # end of kill particle block
            else:
                ip += 1
                ip_m1 += 1

            if ip > np:
                break  # exit loop

        # loops until either np is decreased to ip, or ip increased to np
        if np == 0:
            #  we need at least one particle on the stack
            #  so that the transport routines can exit properly
            stack.np = 1
            e[np_m1] = 0
            iq[np_m1] = 0
            wt[np_m1] = 0


