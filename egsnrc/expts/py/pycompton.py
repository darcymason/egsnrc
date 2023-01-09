from math import log, exp, sqrt
# from egsnrc import random
import random
from egsnrc.constants import RM


def init_random(seed):
    return np.random.default_rng(seed)

def floats_0_1(rng, num):



def py_compton(key, energy):
    """   Subroutine for sampling incoherent (Compton) scattering        """
    # XXX do not use bound, always K-N
    # XXX no e- created yet

    ko = energy / RM # Gamma energy in units of electron rest energy"
    broi = 1 + 2 * ko  # Needed for scattering angle sampling"

    first_time = False
    while True:
        # ... binding effects possible, removed here
        if ko > 2:  # "At high energies the original EGS4 method is most efficient"
            if first_time:
                broi2 = broi * broi
                alph1 = log(broi)
                bro   = 1 / broi
                # alph2 = ko * (broi+1) * bro * bro # not used again, just inline below
                alpha = alph1 + ko * (broi+1) * bro * bro

            while True:  # rejection sampling loop
                key, (rnno15, rnno16, rnno19) = floats_0_1(key, 3)
                if rnno15 * alpha < alph1:  # Use 1/br part
                    br = exp(alph1 * rnno16) * bro
                else:  # Use the br part
                    br = sqrt(rnno16 * broi2 + (1 - rnno16))  * bro

                temp = (1 - br) / (ko * br)
                sinthe = max(0., temp*(2-temp))
                aux = 1 + br * br
                rejf3 = aux - br * sinthe
                if rnno19 * aux <= rejf3:
                    break
        else:  # At low energies it is faster to sample br uniformely
            if first_time:
                bro = 1. / broi
                bro1 = 1 - bro
                rejmax = broi + bro

            while True:
                key, (rnno15, rnno16) = random.floats_0_1(key, 2)
                br = bro + bro1 * rnno15
                temp = (1 - br) / (ko * br)
                sinthe = max(0., temp*(2-temp))
                rejf3 = 1 + br * br - br * sinthe
                if rnno16 * br * rejmax <= rejf3:
                    break
        ]
        first_time = False


        if (bro < br < 1):
            break  [
            # IF( br < 0.99999/broi | br > 1.00001 )
            #     $egs_warning(*,' sampled br outside of allowed range! ',ko,1./broi,br)

    # $RADC_REJECTION

    costhe = 1 - temp
    IF( ibcmp(irl) = 0 ) [ "User wants to use Klein-Nishina, so we are done"
        Uj = 0
        goto :FINISHED-COMPTON-SAMPLING:
    ]

    # FINISHED-COMPTON-SAMPLING
    pesg = br * peig
    # pese = peig - pesg - Uj + prm
    sinthe = sqrt(sinthe)
    # uphi(2,1)
    energy *= br
    # aux = 1 + br*br - 2*br*costhe
    # IF( aux > 1e-8 ) [
    #     costhe = (1-br*costhe)/Sqrt(aux)
    #     sinthe = (1-costhe)*(1+costhe)
    #     IF( sinthe > 0 ) [ sinthe = -Sqrt(sinthe); ]
    #     ELSE [ sinthe = 0; ]
    # ] ELSE [ costhe = 0; sinthe = -1; ]
    # np = np + 1
    # $CHECK-STACK(np,'COMPT')
    # call uphi(3,2)
    # e(np) = pese; iq(np) = -1

    # IF( ibcmpl = 1 | ibcmpl = 3 ) [

    #     " Shell vacancy "
    #     IF( Uj > 1e-3 ) [
    #         edep = pzero

    #         call relax(Uj,shn_array(j),iz_array(j))
    #         "relax will put all particles with energies above ecut,pcut on the "
    #         "stack, the remaining energy will be scored in edep and deposited  "
    #         "locally (via the call to ausgab below)                            "
    #     ]
    #     ELSE [
    #         edep = Uj
    #         edep_local = edep
    #         $AUSCALL($SPHOTONA)
    #     ]
    #     IF( edep > 0 ) [ $AUSCALL($PHOTXAUS); "generates IARG = 4 call" ]

    # ]

    # " Now play Russian Roulette with resulting electrons if the user asked for it"
    # $PLAY RUSSIAN ROULETTE WITH ELECTRONS FROM NPold+1; "TO NP"

    # return

    # :INTERACTION-REJECTED:
    # " Create here a zero energy electron if required (check user codes) "
    # return