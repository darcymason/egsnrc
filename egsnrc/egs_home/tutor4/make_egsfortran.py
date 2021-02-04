#!/usr/bin/env python
import os
import sys
from pathlib import Path
import subprocess
import shutil
import contextlib

# Add a config var for PY_EGS, for modified EGSnrc mortan codes.
# You will also have to modify the .make file for modified
# versions of mortran file stored in this repo (point to PY_EGS location).
HERE  = Path(__file__).resolve().parent
USER_CODE = "tutor4"
CONF = "_linux"
USER_CODE_FORTRAN = f"{USER_CODE}{CONF}.f"
HEN_HOUSE = Path(os.environ["HEN_HOUSE"])
COMPILE_USER_CODE = HEN_HOUSE / "scripts" / "compile_user_code"
PY_EGS = str(HERE.parent.parent / "HEN_HOUSE" / "src") + "/"
LIB_NAME = "egsfortran"
F2PY_OPTIONS = "--quiet --debug"

os.environ['PY_EGS'] = PY_EGS
print(f"Set environment variable PY_EGS to {PY_EGS}")


@contextlib.contextmanager
def working_dir(path):
    save_dir = Path.cwd()
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(save_dir)


with working_dir(HERE):
    print(f"Working directory: {HERE}")
    print(f"Mortran compiling: 'm {USER_CODE}' to {USER_CODE_FORTRAN}")
    command = str(COMPILE_USER_CODE)
    proc = subprocess.run([command, "m", USER_CODE], capture_output=True)
    if proc.returncode != 0:
        print("Error in mortran compile")
        print(proc.stderr)
        sys.exit(-1)

    # Assuming here that 'm' may have done nothing if file existed
    # needs_modification = False
    # with open(USER_CODE_FORTRAN) as f:
    #     if 'subroutine' not in f.readline():
    #         needs_modification = True
    #         f.seek(0)
    #         code = f.read()

    # if needs_modification:
    #     print(f'Modifying {USER_CODE_FORTRAN}  to take out "main"')
    #     with open(USER_CODE_FORTRAN, 'w') as f_out:
    #         f_out.write('      subroutine testxxx\n' + code)


    print(f"Running f2py on {USER_CODE_FORTRAN}")
    proc = subprocess.run(
        f"python3.9 -m numpy.f2py {F2PY_OPTIONS} -c {USER_CODE_FORTRAN} -m {LIB_NAME}".split(),
        # capture_output=True, encoding="utf8"
    )
    if proc.returncode != 0:
        print("Error in numpy.f2py")
        print(f"{LIB_NAME} not created.  Stopping.")
        print(proc.stderr)
    else:
        filenames = HERE.glob(f"{LIB_NAME}*")
        for filename in filenames:
            print(f"Copying {filename.name} shared lib to egsnrc package")
            shutil.copy(filename, HERE.parent.parent)
    # python3.9 -c ";f=
    # python3.9 -m numpy.f2py -c mod_tutor4.f -m egsfortran
    # cp egsfortran* ../../
