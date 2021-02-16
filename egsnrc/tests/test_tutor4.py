import pytest
from pathlib import Path

pytest.importorskip("egsnrc.egsfortran")  # from numpy.f2py, used while in transition
from egsnrc import egsfortran
from egsnrc.egs_home.tutor4 import tutor4
from egsnrc import calcfuncs

import logging
logger = logging.getLogger("egsnrc")



HERE  = Path(__file__).resolve().parent
TEST_DATA = HERE / "data"
TUTOR4_PATH = HERE.parent / "egs_home" / "tutor4" / "tutor4.py"

# @pytest.mark.skipif(sys.platform=="win32")

def known_in_out(filepath, in_types, out_types):
    """Iterator over a filename, yielding known inputs and result"""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    gen = iter(lines)
    for line in gen:
        if not line.startswith("in "):
            continue
        inputs = line[3:].split()  # split after 'in '

        assert len(inputs) == len(in_types)
        inputs = [typ(x.strip()) for x, typ in zip(inputs, in_types)]

        out_line = next(gen)
        while not out_line.startswith("out "):
            out_line = next(gen)

        outputs = out_line[4:].split() # split after 'out '
        if not isinstance(out_types, (list, tuple)):
            out_types = (out_types,)

        assert len(outputs) == len(out_types)
        outputs = [typ(x.strip()) for x, typ in zip(outputs, out_types)]
        if len(outputs) == 1: # original tests only had one output
            yield inputs, outputs[0]
        else:
            yield inputs, outputs


def gen_data_lines(lines):
    """Yield 'data' lines from an egsnrc run"""
    for line in lines:
        if line.find(":") == 35:
            yield line


def line_data(line):
    """Given a `watch` line, return a list of the numbers"""
    data = line.split(":", maxsplit=1)[1]
    data = [float(d.strip()) for d in data.split()]
    return data


def lines_approx_equal(line1, line2, epsilon=0.000002):
    line1 = line1.strip()
    line2 = line2.strip()
    if line1 == line2:
        return True
    else:
        # Allow for small differences in float numbers:
        line1_data = line_data(line1)
        line2_data = line_data(line2)
        if len(line1_data) != len(line2_data):
            return False  # however, should never be the case
        return all(
            abs(d1 - d2) < epsilon
            for d1, d2 in zip(line1_data, line2_data)
        )


class TestTutor4:
    def test_output(self, caplog):
        """Test that Python tutor4 produces known output"""

        logger.propagate = True  # needed for pytest to capture
        caplog.set_level(logging.DEBUG)
        # Ensure proper random initial state
        # (other tests use ranlux)

        egsfortran.init_ranlux(1,0)
        tutor4.main(iwatch=2, high_prec=True)

        std_filename = TEST_DATA / "orig-tutor4-watch2-extra-prec.txt"
        expected = open(std_filename, "r").readlines()
        got = [rec.message.strip('\n') for rec in caplog.records]

        # Test each line of "data" - ignore headings, etc
        iter_got = gen_data_lines(got)
        icount = 0
        for line_expect in gen_data_lines(expected):
            line_got = next(iter_got)
            assert lines_approx_equal(line_expect, line_got)
            icount += 1
        assert icount > 200  # check that test is actually testing lots of lines

    def test_compute_drange(self):
        "Calculate correct values for $COMPUTE-DRANGE in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        # tutor4.init()  # get all data loaded
        # Known inputs for compute-drange from Mortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "compute-drange.txt",
            (int, int, float, float, int, float, float), float
        ):
            # compute_drange(lelec, medium, eke1, eke2, lelke1, elke1, elke2)
            got = calcfuncs.compute_drange(*inputs)
            assert got == pytest.approx(expected,abs=0.0000001)

    def test_calc_tstep(self):
        "Calc correct values for modified $CALCULATE-TSTEP-FROM-DEMFP in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        # tutor4.init()  # get all data loaded
        # Known inputs from Mortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "calc-tstep.txt",
            (int, int, int, int, float, float, float, float, float), float
        ):
            #
            # print("in ", ",".join(str(x) for x in inputs))

            got = calcfuncs.calc_tstep_from_demfp(*inputs)
            # print(got, expected)
            assert got == pytest.approx(expected,abs=0.0000001)

    def test_compute_eloss(self):
        "Calc correct values for $COMPUTE-ELOSS in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        # tutor4.init()  # get all data loaded
        # Known inputs from Mortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "compute-eloss.txt",
            (int, int, float, float, float, int), float
        ):
            #
            # print("in ", ",".join(str(x) for x in inputs))
            got = calcfuncs.compute_eloss(*inputs)

            assert got == pytest.approx(expected,abs=0.0000001)

    def test_compute_eloss_g(self):
        "Calc correct values for $COMPUTE-ELOSS-G in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        # tutor4.init()  # get all data loaded
        # Known input and output from Mortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "compute-eloss-g.txt",
            # lelec, medium, step, eke, elke, lelke, range_
            (int, int, float, float, float, int, float), float
        ):
            #
            # print("in ", ",".join(str(x) for x in inputs))
            got = calcfuncs.compute_eloss_g(*inputs)

            assert got == pytest.approx(expected,abs=0.0000001)

    def test_calculate_xi(self):
        "Calc correct values for $CALCULATE-XI in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        # tutor4.init()  # get all data loaded

        # Need setting here to get to IF conditions where this code applies
        from egsnrc.commons import et_control
        et_control.exact_bca = True

        # Known input and output from Mortran/Fortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "calc-xi.txt",
            # lelec, medium, ekems, rmt2, rmsq, xccl, blccl, step
            (int, int, float, float, float, float, float, float),
            (float, float)
        ):
            #
            # print("in ", ",".join(str(x) for x in inputs))

            got = calcfuncs.calculate_xi(*inputs)
            for a_got, a_expected in zip(got, expected):
                assert a_got == pytest.approx(a_expected,abs=0.0000001)

        et_control.exact_bca = False

    def test_pi_zero(self):
        with pytest.raises(NotImplementedError):
            tutor4.shower(2, 100, 0, 0, 0, 0, 0, 1, 1, 1)