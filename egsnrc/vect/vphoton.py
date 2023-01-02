import numpy


def photon(fparticles, iparticles, howfar, ausgab):
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

    while True:
        # if eig <= pcut[irl_m1]: particle_outcome = PCUT_DISCARD
        iparticles[STATUS,fparticles[ENERGY] < pcut[iparticles[REGION]]] = PCUT_DISCARD

        # if wt == 0.0:  particle_outcome = USER_PHOTON_DISCARD
        wt0_mask = fparticles[WT] == 0
        iparticles[STATUS, wt0_mask] = USER_PHOTON_DISCARD

        alive_mask = (iparticles[STATUS] == 0)  # could use ialive[STATUS] if don't add particles
        falive = fparticles[:, alive_mask]
        ialive = iparticles[:, alive_mask]

        if falive.size == 0:
            break  # all particles have an outcome, we're done this set

        gle = numpy.log(falive[ENERGY])

        #    here to sample no. mfp to transport before interacting

        rnno35 = rng.random(falive.shape[1])
        rnno35[rnno35==0] = 1.0e-30  # could use `numpy.clip` but not quite same as old code
        dpmfp = -numpy.log(rnno35)

        # if medium != 0:
        non_vac_mask = (ialive[MEDIUM] != 0)
        non_vac_media = ialive[MEDIUM, non_vac_mask]

        # set interval gle, ge;
        # lgle is integer index into upcoming arrays
        lgle = (
            ge1[non_vac_media]*gle + ge0[non_vac_media]
        ).astype(numpy.int32)

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

    NUM_PARTICLES = 5

    # float and int state of particles
    # Different energies, one with weight 0
    fphotons = numpy.array(
        [
            (20, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0),
            ( 2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0),
            (18, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0),  # wt= 0
            (12, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0),  # v=w = 1
        ],
        dtype=numpy.float64
    ).transpose()
    # photons in several regions, different mediums
    # could look up medium by region, but easier to just track it
    iphotons = numpy.array(
        [
            (0, 2, 1),
            (0, 1, 1),
            (0, 3, 2),
            (0, 4, 3),
        ],
        dtype=numpy.int32
    ).transpose()

    # Completely made-up data arrays for toy example
    pcut = numpy.array((99, 20, 15, 12, 2), dtype=numpy.float64) # by region
    ge1 = numpy.array((99, 1.1, 2.2, 3.3), dtype=numpy.float64) # by medium
    ge0 = numpy.array((99, 0.1, 0.2, 0.3), dtype=numpy.float64) # by medium

    print("\n\nStarting particles states")
    print(fphotons)
    print(iphotons)

    print("\nStarting 'physics data'")
    print("pcut")
    print(pcut)
    print("ge1, ge0")
    print(ge1)
    print(ge0)

    rng = numpy.random.default_rng(12345)

    photon(fphotons, iphotons, howfar, ausgab)