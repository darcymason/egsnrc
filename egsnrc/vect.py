# vect.py
"""Basic support for vectorized arrays and functions"""

# Particle indices
# ----------------
# For fparticles array (particle float parameters)
PARTICLE_FLOAT_PARAMS = 11
ENERGY, X, Y, Z, U, V, W, WT, DNEAR, TSTEP, USTEP = range(PARTICLE_FLOAT_PARAMS)

# For iparticles array (particle int parameters)
PARTICLE_INT_PARAMS = 3
STATUS, REGION, MEDIUM = range(PARTICLE_INT_PARAMS)
