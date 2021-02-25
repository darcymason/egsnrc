import numpy

MXSTACK = 40  # must be same as when f2py egsfortan created

EPSEMFP: numpy.float64 = 1.E-8
RANDOMIZE_TUSTEP: bool = False
PRESTA_II: numpy.int32 = 0
EM_MACROS_ACTIVE: bool = False
TRANAUSB: numpy.int32 = 0
TRANAUSA: numpy.int32 = 5
MOLLAUSB: numpy.int32 = 8
MOLLAUSA: numpy.int32 = 9
BHABAUSB: numpy.int32 = 10
BHABAUSA: numpy.int32 = 11
ANNIHFAUSB: numpy.int32 = 12
ANNIHFAUSA: numpy.int32 = 13
UPHIAUSB: numpy.int32 = 21
UPHIAUSA: numpy.int32 = 22
BREMAUSB: numpy.int32 = 6
BREMAUSA: numpy.int32 = 7
EGSCUTAUS: numpy.int32 = 1
PEGSCUTAUS: numpy.int32 = 2
ANNIHRAUSB: numpy.int32 = 28
ANNIHRAUSA: numpy.int32 = 14
USERDAUS: numpy.int32 = 3
EIIB: numpy.int32 = 31  # Before EII
EIIA: numpy.int32 = 32  # After EII


EM_MACROS_ACTIVE = False

# XXX extra AUSFL flags for debugging - DLM 2021-02
RANDOMNUM = 40
PRESTAIIA = 41
PRESTAIA = 42
TUSTEPB = 43
HOWFARB = 44
HOWFARA = 45
XYZAUSB = 46
XYZAUSA = 47
UVWAUSB = 48
UVWAUSA = 49

# for EII-Data
MAX_EII_SHELLS = 40  # Maximum number of shells participating
                     # in EII in a simulation
N_EII_BINS = 250  # Number of bins for EII x-section interpolations
MAX_EII_BINS = N_EII_BINS * MAX_EII_SHELLS

# Common block EDGE
MXELEMENT = 100  #  Number of elements
MXSHXSEC = 30  #  Number of shells available
MXSHELL = 6  #  Number of shells treated
MXINTER = 5  #  $MXSHELL-1
MXTRANS = 39  #  Number of possible transitions
MXEDGE = 16  #  max. number of edges above 1 keV
# PHOTOUNIT = i_photo_relax  #  unit number for photo_relax.data
# PHOCSUNIT = i_photo_cs  #  unit number for photo_cs.data