# util.py

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

vals = [1.2345, 22.22, -1e-323]
print(fort_hex(vals))