import numpy as np

REAL = np.float64
ENERGY_PRECISION = np.float32
INTEGER = np.int32
LOGICAL = bool

ARRAY = np.array

test_precision = True  # XXX later should default to False
"""If ``True``, set float constants to match the single precision in default
EGSnrc mortran code.
"""