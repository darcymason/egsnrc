from tqdm import tqdm
import numpy as np
import jax.numpy as jnp


def photon(fparticles, iparticles, howfar, ausgab, max_loops=int(1e3)):
    """Called from `shower` to track photons

    Parameters
    ----------
    fparticles: numpy.ndarray
        2D array of photon states (float64) x NUM_PARTICLES

    iparticles: numpy.ndarray
        2D array of photon states (int32) x NUM_PARTICLES

    howfar:  function
        Callback to `howfar` user-code method

    ausgab:  function
        Callback to `ausgab` user-code method for tracking events and energy

    Returns
    -------
    ircode: int
        1 for a normal return
    """

    iparticles[STATUS] = 0 # set up normal return

    PCUT_DISCARD = 1
    USER_PHOTON_DISCARD = 2
    PAIR_ELECTRONS_KILLED = 3

    # Keep track of particles with no outcome yet, start with all of them
    falive = fparticles
    ialive = iparticles

    for _ in tqdm(range(max_loops)):
        # if eig <= pcut[irl_m1]: particle_outcome = PCUT_DISCARD
        particle_cut_mask = fparticles[ENERGY] < pcut[iparticles[REGION]]
        iparticles[STATUS, particle_cut_mask] = PCUT_DISCARD

        # if wt == 0.0:  particle_outcome = USER_PHOTON_DISCARD
        iparticles[STATUS, fparticles[WT] == 0] = USER_PHOTON_DISCARD

        alive = (iparticles[STATUS] == 0)  # could use ialive[STATUS] if don't add particles
        falive = fparticles[:, alive]
        ialive = iparticles[:, alive]

        if falive.size == 0:
            break  # all particles have an outcome, we're done this set

        gle = np.log(falive[ENERGY])

        #    here to sample no. mfp to transport before interacting

        rnno35 = rng.random(falive.shape[1])
        rnno35[rnno35==0] = 1.0e-30  # could use `numpy.clip` but not quite same as old code
        dpmfp = -np.log(rnno35)

        # if medium != 0:
        non_vac = (ialive[MEDIUM] != 0)
        non_vac_mediums = ialive[MEDIUM, non_vac]
        # set interval gle, ge;
        # lgle is integer index into upcoming arrays
        lgle = (
            ge1[non_vac_mediums]*gle + ge0[non_vac_mediums]
        ).astype(np.int32)

        # xxxx
        # xxxx
        # ******* continue ... ************

        # evaluate gmfpr0 using gmfp(gle)
        # ...
        # ...


def howfar(fparticles, iparticles):
    logger.debug("In howfar")

def ausgab(fparticles, iparticles):
    logger.debug("In ausgab")

if __name__ == "__main__":
    # Toy examples of arrays

    # Particle indexes
    PARTICLE_FLOAT_PARAMS = 11
    ENERGY, X, Y, Z, U, V, W, WT, DNEAR, TSTEP, USTEP = range(PARTICLE_FLOAT_PARAMS)

    PARTICLE_INT_PARAMS = 3
    STATUS, REGION, MEDIUM = range(PARTICLE_INT_PARAMS)

    NUM_PARTICLES = 1e6

    # float and int state of particles
    # Different energies, one with weight 0
    fphotons = np.array(
        [
            (20, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0),
            ( 2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0),
            (18, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),  # wt= 0
            (12, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),  # v=w = 1
        ],
        dtype=np.float64
    ).transpose()
    # photons in several regions, different mediums
    # could look up medium by region, but easier to just track it
    iphotons = np.array(
        [
            (0, 2, 1),
            (0, 1, 1),
            (0, 3, 2),
            (0, 4, 3),
        ],
        dtype=np.int32
    ).transpose()

    num_template_particles = fphotons.shape[1]
    assert iphotons.shape[1] == num_template_particles

    num_repeats_needed = NUM_PARTICLES // num_template_particles

    fphotons = np.repeat(fphotons, num_repeats_needed, axis=1)
    iphotons = np.repeat(iphotons, num_repeats_needed, axis=1)

    # Completely made-up data arrays for toy example
    pcut = np.array((99, 20, 15, 12, 2), dtype=np.float64) # by region
    ge1 = np.array((99, 1.1, 2.2, 3.3), dtype=np.float64) # by medium
    ge0 = np.array((99, 0.1, 0.2, 0.3), dtype=np.float64) # by medium

    print("\n\nStarting particles states")
    print(fphotons)
    print(iphotons)

    print("\nStarting 'physics data'")
    print("pcut")
    print(pcut)
    print("ge1, ge0")
    print(ge1)
    print(ge0)

    rng = np.random.default_rng(12345)

    photon(fphotons, iphotons, howfar, ausgab)