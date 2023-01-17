from math import log, exp, sqrt, sin, cos
from egsnrc import random
from egsnrc.constants import RM, TWO_PI


def py_compton(rng, energy, calc_azimuth=True):
    """Subroutine for sampling incoherent (Compton) scattering"""
    # calc_azimith False used in test_compton only
    # TODO: do not use bound, always K-N
    # TODO: no e- created yet

    ko = energy / RM  # Gamma energy in units of electron rest energy
    broi = 1 + 2 * ko  # Needed for scattering angle sampling

    if ko > 2:  # At high energies the original EGS4 method is most efficient
        broi2 = broi * broi
        alph1 = log(broi)
        bro   = 1 / broi
        # alph2 = ko * (broi+1) * bro * bro # not used again, just inline below
        alpha = alph1 + ko * (broi+1) * bro * bro

        # Set up fake True for first pass through loop
        rnno19 = aux = br = 1; rejf3 = 0.0
        while rnno19 * aux > rejf3 or not (bro < br < 1):  # rejection sampling loop
            _, (rnno15, rnno16, rnno19) = random.floats_0_1(rng, 3)
            if rnno15 * alpha < alph1:  # Use 1/br part
                br = exp(alph1 * rnno16) * bro
            else:  # Use the br part
                br = sqrt(rnno16 * broi2 + (1 - rnno16))  * bro

            temp = (1 - br) / (ko * br)
            sinthe = max(0., temp*(2-temp))
            aux = 1 + br * br
            rejf3 = aux - br * sinthe

            # IF( br < 0.99999/broi | br > 1.00001 )
            #     $egs_warning(*,' sampled br outside of allowed range! ',ko,1./broi,br)

    else:  # At low energies it is faster to sample br uniformly
        bro = 1. / broi
        bro1 = 1 - bro
        rejmax = broi + bro

        # Set up fake True for first pass through loop
        rnno16 = br = 1.0; rejf3 = 0.0
        while rnno16 * br * rejmax > rejf3 or not (bro < br < 1):
            _, (rnno15, rnno16) = random.floats_0_1(rng, 2)
            br = bro + bro1 * rnno15
            temp = (1 - br) / (ko * br)
            sinthe = max(0., temp*(2-temp))
            rejf3 = 1 + br * br - br * sinthe
            # IF( br < 0.99999/broi | br > 1.00001 )
            #     $egs_warning(*,' sampled br outside of allowed range! ',ko,1./broi,br)

    # $RADC_REJECTION

    costhe = 1 - temp
    sinthe = sqrt(sinthe)
    energy *= br
    # Random sample the azimuth
    if calc_azimuth:
        _, (azimuth_ran,) = random.floats_0_1(rng, 1)
        phi = TWO_PI * azimuth_ran
        sinphi = sin(phi)
        cosphi = cos(phi)
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

    return energy, sinthe, costhe, sinphi, cosphi


if __name__ == "__main__":
    random.set_array_library("numpy")
    rng = random.initialize(42)
    import numpy as np

    for in_energy in (0.01, 0.1, 1, 10):  # 0.1, 1, 10
        print("-----------------------------------------------------------")
        print(f"Sample outputs for energy {in_energy} MeV:")
        print("\t".join("energy sinthe costhe sinphi cosphi sin2+cos2the check_E'".split()))
        for i in range(10):
            result = py_compton(rng, in_energy)
            print("\t".join(f"{result[i]:0.4f}" for i in range(len(result))), end="")
            sinthe, costhe, sinphi, cosphi = result[1:]
            checkthe = sinthe**2 + costhe**2
            check_e = in_energy / (1 + in_energy/RM*(1-costhe))
            print("\t" + "\t".join((f"{checkthe:.6f}", f"{check_e:.6f}")))

    # print("Sample angular results to compare with known distribution")

    # from math import acos, degrees
    # energies = np.array((0.01, 0.1, 1, 10))
    # histo_deg = np.zeros((180, len(energies)), dtype=np.int32)
    # for i, in_energy in enumerate(energies):
    #     for _ in range(80_000):
    #         energy, _, costhe = py_compton(rng, in_energy)[:3]
    #         theta_bin = int(degrees(acos(costhe)))
    #         histo_deg[theta_bin, i] += 1

    # from contextlib import redirect_stdout
    # with open(r"c:\temp\pycompt.tsv", 'w') as f:
    #     with redirect_stdout(f):
    #         print("\t".join(("Angle",*(str(e) for e in energies))))
    #         for deg in range(180):
    #             counts = [f"{x:5}" for x in histo_deg[deg]]
    #             print("\t".join((f"{deg:3}", *counts)))


    from time import perf_counter

    num_part = 100_000
    print(f"pycompton: starting {num_part} 'particles'")
    energy = 1
    start = perf_counter()
    for i in range(num_part):
        py_compton(rng, energy, calc_azimuth=False)
    end = perf_counter()
    print(f"Done in {(end - start):.5} seconds")