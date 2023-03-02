from egsnrc.config import device_jit
from egsnrc.particles import replace_e_uvw


# XXX photo() is a placeholder for now, obviously trivial for photons,
# but will need for electron later
@device_jit
def photo(gid, rng_states, p):
    """Return a modified particle from photoelectric interaction

    Returns
    -------
    p   Particle
        The modified particle with no energy
    """

    # For now, re-use existing e_uvw replace to just replace energy
    return replace_e_uvw(p, 0.0, p.u, p.v, p.w)
