import numpy as np
from egsnrc import random
from egsnrc.constants import RM, TWO_PI


def pyvect_compton(rng, energy, calc_azimuth=True):
    """Subroutine for sampling Compton scattering using numpy vectors

    Params
    ------
    rng:    random generator
    energy: np.ndarray
        Array of energies to calculate
    calc_azimith:  bool
        Whether to random sample azimuthal angle.
        False for some tests to stay in synch with random sequence from mortran run

    Returns
    -------
    arrays of energy, sinthe, costhe, sinphi, cosphi
    """
    # calc_azimith False used in test_compton only
    # TODO: do not use bound, always K-N
    # TODO: no e- created yet

    ko = energy / RM  # Gamma energy in units of electron rest energy
    broi = 1 + 2 * ko  # Needed for scattering angle sampling
    bro   = 1 / broi

    # Final values passing rejection sampling, for last calcs before return
    # Could instead update energy directly rather than through
    all_sinthe = np.empty_like(energy)
    all_temp = np.empty_like(energy)
    all_br = np.empty_like(energy)

    # ko > 2: At high energies the original EGS4 method is most efficient
    ko_gt2_mask = ko > 2
    i_notdone = ko_gt2_mask.nonzero()[0]  # indices we still need to process
    if len(i_notdone):
        broi2 = np.square(broi, where=ko_gt2_mask)
        alph1 = np.log(broi, where=ko_gt2_mask)
        # alph2 = ko * (broi+1) * bro * bro # not used again, just inline below
        alpha = alph1 + np.where(ko_gt2_mask, ko * (broi + 1) * bro * bro, 0)
        while len(i_notdone):
            _, (rnno15, rnno16, rnno19) = random.floats_0_1(rng, (3, len(i_notdone)))
            alph1_notdone = alph1[i_notdone]
            bro_notdone = bro[i_notdone]
            br_condn = rnno15 * alpha[i_notdone] < alph1_notdone
            # Use either 1/br part or br part depending on condition
            br = np.where(
                br_condn,
                np.exp(alph1_notdone * rnno16, where=br_condn) * bro_notdone, # 1/br part
                np.sqrt(rnno16 * broi2[i_notdone] + (1 - rnno16), where=~br_condn) * bro_notdone
            )
            # rejection sampling
            temp = (1 - br) / (ko[i_notdone] * br)
            sinthe = temp * (2 - temp)
            sinthe.clip(min=0.0, out=sinthe)
            aux = 1 + np.square(br)
            rejf3 = aux - br * sinthe
            pass_rej = rnno19 * aux <= rejf3
            newlydone_mask = pass_rej & (br < 1) & (br > bro_notdone)

            # Calc indices of global arrays to update
            i_newlydone = i_notdone[newlydone_mask]

            # Update the global arrays with elements that passed in this loop
            all_br[i_newlydone] = br[newlydone_mask]
            all_temp[i_newlydone] = temp[newlydone_mask]
            all_sinthe[i_newlydone] = sinthe[newlydone_mask]

            # Update which particle indexes have been completed
            i_notdone = i_notdone[~newlydone_mask]

    # ko <= 2: At low energies it is faster to sample br uniformly
    ko_le2_mask = ~ko_gt2_mask
    i_notdone = ko_le2_mask.nonzero()[0]  # indices still need to process
    if len(i_notdone):
        # Don't use where clause as these are fast
        bro1 = 1 - bro
        rejmax = broi + bro

        # Set all of the current ones as 'not done'
        while len(i_notdone):
            _, (rnno15, rnno16) = random.floats_0_1(rng, (2, len(i_notdone)))
            bro_notdone = bro[i_notdone]
            br = bro_notdone + bro1[i_notdone] * rnno15
            temp = (1 - br) / (ko[i_notdone] * br)
            sinthe = temp * (2 - temp)
            sinthe.clip(min=0.0, out=sinthe)
            rejf3 = 1 + np.square(br) - br * sinthe
            pass_rej = rnno16 * br * rejmax[i_notdone] <= rejf3
            newlydone_mask = pass_rej & (br < 1) & (br > bro_notdone)

            # Calc indices of global arrays to update
            i_newlydone = i_notdone[newlydone_mask]

            # Update the global arrays with elements that passed in this loop
            all_br[i_newlydone] = br[newlydone_mask]
            all_temp[i_newlydone] = temp[newlydone_mask]
            all_sinthe[i_newlydone] = sinthe[newlydone_mask]

            # Update which particle indexes have been completed
            i_notdone = i_notdone[~newlydone_mask]

            # IF( br < 0.99999/broi | br > 1.00001 )
            #     $egs_warning(*,' sampled br outside of allowed range! ',ko,1./broi,br)

    # $RADC_REJECTION

    result_costhe = 1 - all_temp
    result_sinthe = np.sqrt(all_sinthe)
    result_energy = energy * all_br
    # Random sample the azimuth
    if calc_azimuth:
        _, (azimuth_ran,) = random.floats_0_1(rng, len(energy))
        phi = TWO_PI * azimuth_ran
        sinphi = np.sin(phi)
        cosphi = np.cos(phi)
    else:
        sinphi = cosphi = 99
    # aux = 1 + br*br - 2*br*costhe
    # if aux > 1e-8:
    #     costhe = (1-br*costhe)/sqrt(aux)
    #     sinthe = (1-costhe)*(1+costhe)
    #     sinthe = -sqrt(sinthe) if sinthe > 0 else 0
    # else:
    #     costhe = 0
    #     sinthe = -1

    return result_energy, result_sinthe, result_costhe, sinphi, cosphi


if __name__ == "__main__":
    if False:
        random.set_array_library("sequence")
        rand_sequence = ([
            [
                [0.22987532615661621, 0.78682494163513184],
                [0.76323693990707397, 0.19077306985855103]
            ],
            [
                [0.55831658840179443,  4.9074590206146240E-002]
            ],
            [
                [0.55832588672637939, 0.29218196868896484]
            ],
            [
                [0.65860986709594727,  4.0171027183532715E-002]
            ],
            [
                [1.7512619495391846E-002,  0.44601458311080933]
            ],
            [
                [0.91792672872543335, 0.54524618387222290],
                [0.28818345069885254,  3.5068333148956299E-002]
            ],
            [
                [0.13349604606628418, 0.85515218973159790],
                [0.54984831809997559,  5.7695209980010986E-002],
            ],
            [
                [0.49270433187484741, 0.94706720113754272],
                [0.84880810976028442, 0.67912405729293823],
                [0.30815130472183228, 0.37652158737182617],
                [0.96473675966262817, 0.67000657320022583],
                [0.96202033758163452, 0.26576608419418335],
            ],
            [
                [0.62973368167877197, 0.13346099853515625],
            ],
        ])
        energies = np.ones(9)
        result_e_costhe = [  # NOTE this is "costhe before changes" i.e. before UPHI etc.
            (0.81141922919270071, 0.88123946893996008),
            (0.64820104040175874, 0.72266489438445403),
            (0.64820844647776044, 0.72267390145936683),
            (0.72808421050232852, 0.80915849456905908),
            (0.21745297830124174, -0.83892959634582454),
            (0.43304114464512694, 0.33097491664454248),
            (0.64145609910772927, 0.71437552426368434),
            (0.96974936034766024, 0.98405975176254867),
            (0.70508445157959732, 0.78626455389275030),
        ]

        rng = random.initialize(rand_sequence, vect=True)
        result = pyvect_compton(rng, energies, calc_azimuth=False)
        energy, _, costhe = result[:3]

        import pytest
        for e, cos, (expect_e, expect_cos) in zip(energy, costhe, result_e_costhe):
            assert e == expect_e
            assert cos == expect_cos
        print("Assertions passed")
    else:
        from timeit import timeit
        random.set_array_library("numpy")
        rng = random.initialize(42)
        _, energies = random.floats_0_1(rng, 100_000)
        energies = np.array(energies) + 0.001

        from time import perf_counter
        start = perf_counter()
        print(f"PyVECTcompton: Starting {len(energies)} 'particles'")
        pyvect_compton(rng, energies, calc_azimuth=False)
        end = perf_counter()
        print(f"Done in {(end - start):.5} seconds")