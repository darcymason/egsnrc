# #############################################################################
#
#   EGSnrc ranlux random number generator
#   Copyright (C) 2015 National Research Council Canada
#
#   This file is part of EGSnrc.
#
#   EGSnrc is free software: you can redistribute it and/or modify it under
#   the terms of the GNU Affero General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   EGSnrc is distributed in the hope that it will be useful, but WITHOUT ANY
#   WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
#   FOR A PARTICULAR PURPOSE.  See the GNU Affero General Public License for
#   more details.
#
#   You should have received a copy of the GNU Affero General Public License
#   along with EGSnrc. If not, see <http://www.gnu.org/licenses/>.
#
# #############################################################################
#
#   Authors:         Iwan Kawrakow, 2000
#                    Dave Rogers, 2000
#
#   Contributors:
#
# #############################################################################
#
#   Parts of the EGS code that originate from SLAC are distributed by NRC
#   under the terms of the AGPL 3.0 licence, in agreement with SLAC.
#
#   The contributors named above are only those who could be identified from
#   this file's revision history.
#
#   A large number of people have been involved with the development of EGS
#   over the years. Many details are in the manual. The authors want to point
#   out the central role of Ralph Nelson of SLAC in the development of EGS
#   over decades. We all owe Ralph a huge debt of gratitude.  Similarly, the
#   role of Alex Bielajew while he was at NRC was critical in many aspects of
#   code development. Hideo Hirayama was involved with the development of EGS4
#   while at SLAC and with Yosh Namito and Syuichi Ban at KEK has continued
#   developments on EGS4. As well many others and users from around the world
#   have assisted in developing and making available the code, in particular
#   Robert D. Stewart.
#
# #############################################################################
#
#   Subtract-and-borrow random number generator proposed by Marsaglia and
#   Zaman, implemented by F. James with the name RCARRY in 1991, and later
#   improved by Martin Luescher in 1993 to produce 'Luxury Pseudorandom
#   Numbers'. Fortran 77 coded by F. James, 1993.
#
#   References:
#
#   M. Luscher, Computer Physics Communications  79, 100 (1994)
#   F. James, Computer Physics Communications 79, 111 (1994)
#
#   The following 5 'luxury levels' can be used:
#
#   level 0: equivalent to the original RCARRY of Marsaglia and Zaman, very
#            long period, but fails many tests.
#
#   level 1: considerable improvement in quality over level 0, now passes the
#            gap test, but still fails spectral test.
                                                                              "
#   level 2: passes all known tests, but theoretically still defective.
                                                                              "
#   level 3: any theoretically possible correlations have very small chance of
#            being observed.
#
#   level 4: highest possible luxury, all 24 bits chaotic.
#
#   Note that in actual EGSnrc calculations we have obtained incorrect results
#   when using level 0, but level 1 and higher give the same results.
#
#   If the generator is not initialized explicitly, it will be with the
#   default luxury level defined in DEFAULT_LL (which is 1 by default).
#
# #############################################################################

import numpy as np


DEFAULT_LL: np.int32 = 1
not_initialized: bool = True
nskipll = np.array([0,24,73,199,365])
icon = 2147483563

def ranlux(rng_array: np.float64):
    # " Input/print(                  "
    # "==============================="
    # $REAL       rng_array(24)
    # integer*4   seedin,luxury_level
    # integer*4   state(25)
    # $INTEGER    ounit
    # character*(*) fmt_flags

    # " Internal variables            "
    # "==============================="

    # integer*4   seeds(24),carry;    " The state of the generator "
    # integer*4   i24,j24;            " The rng seeds "
    # integer*4   next(24);           " for convinience "

    # integer*4   jseed_dflt,nskip,icon,j,k,status,jseed,nskipll(0:4),icarry
    # logical     not_initialized
    # real*4      twom24,twop24
    # integer*4   uni

    # save        seeds,carry,i24,j24,next,twom24,not_initialized,
    #             nskip,twop24,nskipll


    jseed_dflt/314159265/,


    if not_initialized:
        not_initialized = False
        nskip = nskipll(DEFAULT_LL)
        twom24 = 1
        twop24 = 1
        jseed = jseed_dflt
        for j in range(1, 24+1):
            twom24 = twom24 * 0.5
            twop24 = twop24 * 2
            k = jseed / 53668
            jseed = 40014 * (jseed - k*53668) - k*12211
            if jseed < 0:
                jseed = jseed + icon
            seeds[j] = jseed % 16777216
            next_[j] = j-1
        next_[1] = 24
        i24 = 24
        j24 = 10
        carry = 0
        if seeds[24] == 0:
            carry = 1

    for j in range(24):
        uni = seeds[j24] - seeds[i24] - carry
        if  uni < 0:
            uni = uni + 16777216
            carry = 1
        else:
            carry = 0
        seeds[i24] = uni
        #  if  uni == 0 [ uni = twom24*twom24; ]
        i24 = next_[i24]
        j24 = next_[j24]
        if uni >= 4096:
            rng_array[j] = uni * twom24
        else:
            rng_array[j] = uni * twom24 + seeds[j24]*twom24*twom24

    if nskip > 0:
        for j in range(1,nskip+1):
            uni = seeds[j24] - seeds[i24] - carry
            if  uni < 0:
                uni = uni + 16777216
                carry = 1
            else:
                carry = 0

    seeds[i24] = uni
    i24 = next_[i24]
    j24 = next_[j24]
    return


def init_ranlux(luxury_level, seedin):
    jseed = seedin
    if  jseed <= 0:
        jseed = jseed_dflt
    if  luxury_level < 0 or luxury_level > 4:
        luxury_level = DEFAULT_LL
    nskip = nskipll(luxury_level)

    logger.info(
        ' ***************** RANLUX initialization ******************\n'
        f' luxury level: {luxury_level}\n'
        f' initial seed: {jseed}\n'
        '***********************************************************'
    )

    not_initialized = False
    twom24 = 1; twop24 = 1
    for j in range(1,25) [
        twom24 = twom24 * 0.5
        twop24 = twop24 * 2
        k = jseed / 53668
        jseed = 40014*(jseed - k*53668) - k*12211
        if jseed < 0:
            jseed = jseed + icon
        seeds[j] = jseed % 16777216  # mod()
        next_[j] = j - 1
    ]
    next_[1] = 24
    i24 = 24
    j24 = 10
    carry = 0
    if seeds[24] == 0:
        carry = 1

    return

def get_ranlux_state(state):
    for j in range(1,24+1):
        state[j] = seeds[j]
    state[25] = i24 + 100*(j24 + 100*nskip)
    if  carry > 0:
        state[25] = -state[25]
return

def set_ranlux_state(state):
    twom24 = 1
    twop24 = 1
    for j in range(1,24+1):
        twom24 = twom24 * 0.5
        twop24 = twop24 * 2
        next_[j] = j - 1
    next_[1] = 24
    for j in range(1,24+1):
        seeds[j] = state[j]
    if state[25] <= 0:
        status = -state[25]
        carry = 1
    else:
        status = state[25]
        carry = 0
    nskip = status / 10000
    status = status - nskip*10000
    j24 = status/100
    i24 = status - 100*j24
    if  j24 < 1 or j24 > 24 or i24 < 1 or i24 > 24:
        msg = (
            "**** Error in set_ranlux_state: seeds outside of allowed range!"
            f"status = {state[25]}, nskip = {nskip}, i24 = {i24}, j24 = {j24}"
        raise ValueError(msg)
    not_initialized = False
    return

def show_ranlux_seeds(ounit):

if  carry > 0 ) [ icarry = 1; ]
else:
    icarry = 0; ]
write(ounit,'(a,i4,a,2i3,a,i2,$)')
  ' skip = ',nskip,' ix jx = ',i24,j24,' carry = ',icarry
return

def print_ranlux_seeds(ounit,fmt_flags):

if  carry > 0 ) [ icarry = 1; ]
else:
    icarry = 0; ]
write(ounit,fmt_flags) nskip,i24,j24,icarry
return

end

;" ************************* end of ranlux.mortran  *************************"
