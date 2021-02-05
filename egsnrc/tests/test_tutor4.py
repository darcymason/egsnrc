import pytest
from pathlib import Path

pytest.importorskip("egsnrc.egsfortran")  # from numpy.f2py, used while in transition
from egsnrc import egsfortran
from egsnrc.egs_home.tutor4 import tutor4

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
        if not out_line.startswith("out "):
            raise ValueError("'out' line must follow 'in' line")
        outputs = out_line[4:].split() # split after 'out '
        if not isinstance(out_types, (list, tuple)):
            out_types = (out_types,)

        assert len(outputs) == len(out_types)
        outputs = [typ(x.strip()) for x, typ in zip(outputs, out_types)]
        if len(outputs) == 1: # original tests only had one output
            yield inputs, outputs[0]
        else:
            yield inputs, outputs

class TestTutor4:
    def xxxtest_output(self, capfd):
        """Test that (partially) Python tutor4 produces known output"""

        # Ensure proper random initial state
        # (other tests use ranlux)
        egsfortran.init_ranlux(1,0)
        tutor4.main()
        captured = capfd.readouterr()

        # Test last line of history 10
        lines = captured.out.splitlines()
        rev_lines_iter = reversed(lines)
        while not next(rev_lines_iter).startswith(" END OF HISTORY      10"):
            pass
        last_line_hist10 = next(rev_lines_iter)

        # Could test more, but testing last couple of lines of last history
        # should pretty much guarantee had exact 'trajectory' through
        # random numbers and all physics as original EGSnrc mortran/fortran
        # (within the significant figures shown)
        assert (
            "Discard -user request             :    1   "
            "16.496  -1   3  -0.004  -0.001   0.100 -0.280"
            "  0.259  0.924         0 1.000E+00"
        ) in last_line_hist10
        # Ignore test outputs we might be generating
        while (secondlast := next(rev_lines_iter)).startswith(("in", "out", "fn:")):
            pass
        expected_2nd_last_hist10 = (
            "17.150  -1   2   0.006  -0.005   0.063 -0.240 -0.009  0.971"
        )
        assert expected_2nd_last_hist10 in secondlast

    def test_compute_drange(self):
        "Calculate correct values for $COMPUTE-DRANGE in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        tutor4.init()  # get all data loaded
        # Known inputs for compute-drange from Mortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "compute-drange.txt",
            (int, int, float, float, int, float, float), float
        ):
            # compute_drange(lelec, medium, eke1, eke2, lelke1, elke1, elke2)
            got = tutor4.compute_drange(*inputs)
            assert got == pytest.approx(expected,abs=0.0000001)

    def test_calc_tstep(self):
        "Calc correct values for modified $CALCULATE-TSTEP-FROM-DEMFP in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        tutor4.init()  # get all data loaded
        # Known inputs from Mortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "calc-tstep.txt",
            (int, int, int, int, float, float, float, float, float), float
        ):
            #
            # print("in ", ",".join(str(x) for x in inputs))

            got = tutor4.calc_tstep_from_demfp(*inputs)
            # print(got, expected)
            assert got == pytest.approx(expected,abs=0.0000001)

    def test_compute_eloss(self):
        "Calc correct values for $COMPUTE-ELOSS in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        tutor4.init()  # get all data loaded
        # Known inputs from Mortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "compute-eloss.txt",
            (int, int, float, float, float, int), float
        ):
            #
            # print("in ", ",".join(str(x) for x in inputs))
            got = tutor4.compute_eloss(*inputs)

            assert got == pytest.approx(expected,abs=0.0000001)

    def test_compute_eloss_g(self):
        "Calc correct values for $COMPUTE-ELOSS-G in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        tutor4.init()  # get all data loaded
        # Known input and output from Mortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "compute-eloss-g.txt",
            # lelec, medium, step, eke, elke, lelke, range_
            (int, int, float, float, float, int, float), float
        ):
            #
            # print("in ", ",".join(str(x) for x in inputs))
            got = tutor4.compute_eloss_g(*inputs)

            assert got == pytest.approx(expected,abs=0.0000001)

    def xxxtest_calculate_xi(self):
        "Calc correct values for $CALCULATE-XI in Python"
        # Compare against ones captured from TUTOR4 run with extra prints
        tutor4.init()  # get all data loaded

        # Need setting here to get to IF conditions where this code applies
        from egsnrc.commons import et_control
        et_control.exact_bca = True

        # Known input and output from Mortran/Fortran tutor4 run
        for inputs, expected in known_in_out(TEST_DATA / "calc-xi.txt",
            # medium, ekems, rmt2, rmsq, xccl, blccl, step
            (int, float, float, float, float, float, float), (float, float)
        ):
            #
            # print("in ", ",".join(str(x) for x in inputs))


            got = tutor4.calculate_xi(*inputs)
            if isinstance(got, (list, tuple)):
                for a_got, a_expected in zip(got, expected):
                    assert a_got == pytest.approx(a_expected,abs=0.0000001)
            else:
                assert got == pytest.approx(expected,abs=0.0000001)

        et_control.exact_bca = False