import pytest
from egsnrc import config  # import and line below must precede calcfuncs
config.test_precision = True

from pathlib import Path
pytest.importorskip("egsnrc.egsfortran")  # from numpy.f2py, used while in transition
from egsnrc import egsfortran
from egsnrc.egs_home.tutor4 import tutor4
from egsnrc import calcfuncs

from egsnrc.util import float_from_fort_hex as float_hex

import logging
logger = logging.getLogger("egsnrc")



HERE  = Path(__file__).resolve().parent
TEST_DATA = HERE / "data"
TUTOR4_PATH = HERE.parent / "egs_home" / "tutor4" / "tutor4.py"

# @pytest.mark.skipif(sys.platform=="win32")

def known_in_out(filepath, in_types, out_types, description=""):
    """Iterator over a filename, yielding known inputs and result"""
    with open(filepath, 'r') as f:
        lines = f.readlines()

    gen = iter(lines)
    in_linestart = "in " + description
    out_linestart = "out " + description
    for line in gen:
        if not line.startswith(in_linestart):
            continue
        inputs = line[len(in_linestart):].split()  # split after 'in '

        assert len(inputs) == len(in_types), "Mismatch in input types and inputs"
        inputs = [typ(x.strip()) for x, typ in zip(inputs, in_types)]

        out_line = next(gen)
        while not out_line.startswith(out_linestart):
            out_line = next(gen)

        outputs = out_line[len(out_linestart):].split() # split after 'out '
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


def check_known_in_out(
    filename, func, input_types, output_types, description, min_count
):
    """Compare against ones captured from TUTOR4 run with extra prints
    """
    icount = 0
    for inputs, expected in known_in_out(
        filename, input_types, output_types, description
    ):
        # print("in ", ",".join(str(x) for x in inputs))
        output = func(*inputs)
        if isinstance(expected, (list, tuple)):
            for got, expect in zip(output, expected):
                assert expect == got
        else:
            assert expected == output
        icount += 1
    assert icount > min_count

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
        match = all(
            abs(d1 - d2) < epsilon
            for d1, d2 in zip(line1_data, line2_data)
        )
        return match


class TestTutor4:
    """Tests related to inputs/outputs of EGSnrc Tutor4 example simulation"""
    def setup(self):
        tutor4.init()
        egsfortran.init_ranlux(1,0)
        egsfortran.ranlux(egsfortran.randomm.rng_array)
        egsfortran.randomm.rng_seed = 1

    def test_output_watch2_elec(self, caplog):
        """Test Python tutor4 produces known output - electrons"""
        logger.propagate = True  # needed for pytest to capture
        caplog.set_level(logging.DEBUG)
        # Ensure proper random initial state
        # (other tests use ranlux)

        # Expected data
        std_filename = TEST_DATA / "fort_tut4_elec_outputs.txt"
        with open(std_filename, "r") as f:
            expected = f.readlines()

        tutor4.main(iqin=-1, iwatch=2, high_prec=True, ncase=20)
        got = [rec.message.strip('\n') for rec in caplog.records]

        # Test each line of "data" - ignore headings, etc
        iter_expect = gen_data_lines(expected)
        icount = 0
        for line_got in gen_data_lines(got):
            line_expect = next(iter_expect)
            assert lines_approx_equal(line_expect, line_got), (
                f"Expected line:\n{line_expect}\ngot:\n{line_got}"
            )
            icount += 1
        assert icount > 400  # check that test is actually testing lots of lines

    def test_output_watch2_pos(self, caplog):
        """Test Python tutor4 produces known output - electrons"""
        logger.propagate = True  # needed for pytest to capture
        caplog.set_level(logging.DEBUG)
        # Ensure proper random initial state
        # (other tests use ranlux)

        # Expected data
        std_filename = TEST_DATA / "fort_tut4_pos_outputs.txt"
        with open(std_filename, "r") as f:
            expected = f.readlines()
        tutor4.main(iqin=+1, iwatch=2, high_prec=True, ncase=20)

        got = [rec.message.strip('\n') for rec in caplog.records]

        # Test each line of "data" - ignore headings, etc
        iter_expect = gen_data_lines(expected)
        icount = 0
        for line_got in gen_data_lines(got):
            line_expect = next(iter_expect)
            assert lines_approx_equal(line_expect, line_got), (
                f"Expected line:\n{line_expect}\ngot:\n{line_got}"
            )
            icount += 1
        assert icount > 400  # check that test is actually testing lots of lines

    def test_compute_drange(self):
        "Calculate correct values for $COMPUTE-DRANGE in Python"

        check_known_in_out(
            TEST_DATA / "fort_tut4_calcfuncs.txt",
            calcfuncs.compute_drange,
            (int, int, float_hex, float_hex, int, float_hex, float_hex),
            float_hex,
            "compute-drange:",
            600
        )

    def test_calc_tstep(self):
        "Calc correct values for modified $CALCULATE-TSTEP-FROM-DEMFP in Python"
        # tutor4.init()  # need arrays etc. loaded, so if this is on its own
        # Somehow vacdst was changed while running other tests, to
        #   just under 1e8 (99_999_999.9xxx), so reset for this test
        #   (used in one branch of calc_tstep_from_demfp)
        egsfortran.bounds.vacdst = 1e8
        check_known_in_out(
            TEST_DATA / "fort_tut4_calcfuncs.txt",
            calcfuncs.calc_tstep_from_demfp,
            (int,)*4 + (float_hex,)*5,
            float_hex,
            "calc-tstep-from-demfp:",
            75
        )

    def test_compute_eloss(self):
        "Calc correct values for $COMPUTE-ELOSS in Python"
        check_known_in_out(
            TEST_DATA / "fort_tut4_calcfuncs.txt",
            calcfuncs.compute_eloss,
            (int, int, float_hex, float_hex, float_hex, int),
            float_hex,
            "compute-eloss:",
            500
        )

    def test_compute_eloss_g(self):
        "Calc correct values for $COMPUTE-ELOSS-G in Python"
        check_known_in_out(
            TEST_DATA / "fort_tut4_calcfuncs.txt",
            calcfuncs.compute_eloss_g,
            (int, int, float_hex, float_hex, float_hex, int, float_hex),
            float_hex,
            "compute-eloss-g:",
            500
        )

    def test_calculate_xi(self):
        "Calc correct values for $CALCULATE-XI in Python"

        check_known_in_out(
            TEST_DATA / "calc_not_exact_bca.txt",
            calcfuncs.calculate_xi,
            (int, int) + (float_hex,)*6,
            (float_hex, float_hex),
            "calc-xi:",
            80
        )

    def test_pi_zero(self):
        with pytest.raises(NotImplementedError):
            callbacks = {x: None for x in ('hownear', 'howfar', 'ausgab')}
            tutor4.shower(2, 100, 0, 0, 0, 0, 0, 1, 1, 1, callbacks)
