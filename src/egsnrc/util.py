# util.py
from numba import cuda


def cuda_details():
    try:
        from cuda.cuda import (
            CUdevice_attribute, cuDeviceGetAttribute, cuDeviceGetName, cuInit
        )
    except ImportError as e:
        return (
            f"{str(e)}\n"
            "** GPU details not available\n"
            "In Colab, use `!pip install cuda-python` to see GPU specs\n"
        )

    # Initialize CUDA Driver API
    (err,) = cuInit(0)

    # Get attributes
    err, DEVICE_NAME = cuDeviceGetName(128, 0)
    DEVICE_NAME = DEVICE_NAME.decode("ascii").replace("\x00", "")

    attrs = {'DEVICE_NAME': DEVICE_NAME.strip()}
    attr_names = "MAX_THREADS_PER_BLOCK MAX_GRID_DIM_X MULTIPROCESSOR_COUNT".split()
    for attr in attr_names:
        err, attrs[attr] =  cuDeviceGetAttribute(
        getattr(CUdevice_attribute, f"CU_DEVICE_ATTRIBUTE_{attr}"), 0
    )

    return attrs


def float_from_fort_hex(hex_str):
    return float.fromhex(py_hex_from_fort_hex(hex_str))

def py_hex_from_fort_hex(hex_str) -> str:
    """Return a Python float from a Fortran Z17 exact float hex output"""
    # See e.g. https://en.wikipedia.org/wiki/Double-precision_floating-point_format
    # for description.  52 bit fraction is stored, but could be implied leading
    # zero or 1.
    hex_str = hex_str.strip()
    if hex_str == "0":
        return '0x0.0p+0'
    if len(hex_str) < 16:
        hex_str = f"{hex_str:>016}"

    # first three hex chars = 1 bit sign, 11-bit exponent
    sign_exp, fraction = int(hex_str[:3], 16), int(hex_str[3:], 16)
    sign = sign_exp >> 11
    sign_char = "-" if sign == 1 else " "

    exponent = sign_exp & 0x7FF
    leading_bit = "1"
    if exponent == 0:
        leading_bit = "0"
        exponent = 1
    fraction_hex = hex(fraction)[2:]  # exclude the leading "0x" Python adds
    return f"{sign_char}0x{leading_bit}.{fraction_hex:>013}p{exponent-1023}"

def fort_hex(values):
    """Return a list of floats in Fortran hex format

    e.g. 2 -> 4000000000000000;  -2 -> C000000000000000

    """

    if not isinstance(values, (list, tuple)):
        values = [values]

    # Python float.hex gives values like 0x1.3C083126E978Dp+0 for 1.2345
    # toss the 0x1, convert part after p to +1023, then flip high bit if negative
    py_h = [float(num).hex().upper() for num in values]
    parts = [
        (0x800 if h.startswith("-") else 0, h.find("X"), h.find("P"))
        for h in py_h
    ]
    fort_h = [
        f"{hex(int(h[p+1:])+1023+sgn)[2:].upper()}{h[x+3:p]}"
        for h, (sgn, x, p) in zip(py_h, parts)
    ]
    # Change "3FF" to "        0" to match Fortran
    fort_h = [
        f"{'0':>16}" if h.endswith("3FF0") else h
        for h in fort_h
    ]
    return " " + " ".join(fort_h)


def for_E18(values):
    """High precision exp format which give 18 sig figs for any number

    This is a debugging aid for comparing Python values to Fortran code
    """
    if not isinstance(values, (list, tuple)):
        values = [values]

    py_e = [f"{num:.17E}" for num in values]  # e.g. 1.23456...E-04
    parts = [
        (e[0] if e.startswith("-") else " ", e.find("."), e.find("E"))
        for e in py_e
    ]
    for_e = [
        f"{sgn}0.{e[dot-1]}{e[dot+1:iE]}E{int(e[iE+1:])+1:+03d}"
        for e, (sgn, dot, iE) in zip(py_e, parts)
    ]

    return " " + " ".join(for_e)

# vals = [1.3E-4, -1234567890.123456789, -23.4, 0.0000000123]
# print([f"{x:.18e}" for x in vals])
# print(for_E18(vals))


# from https://towardsdatascience.com/cuda-by-numba-examples-7652412af1ee
# Use like:
# with CUDATimer(stream) as cudatimer:
#    kernel_func[blocks_per_grid, threads_per_block, stream](dev_a, dev_a_reduce)
# print(f"Elapsed time {cudatimer.elapsed:.2f} ms")
class CUDATimer:
    def __init__(self, stream=0):
        self.stream = stream
        self.event = None  # in ms

    def __enter__(self):
        self.event_beg = cuda.event()
        self.event_end = cuda.event()
        self.event_beg.record(stream=self.stream)
        return self

    def __exit__(self, type, value, traceback):
        self.event_end.record(stream=self.stream)
        self.event_end.wait(stream=self.stream)
        self.event_end.synchronize()
        self.elapsed = self.event_beg.elapsed_time(self.event_end)
