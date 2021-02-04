.. _tutor:

Tutor Examples
==============

The tutor4 example from EGSnrc is being used as a testbed for 
checking ongoing switching of code to Python.  The Numpy ``f2py``
module is used to make the original compiled Fortran code (for tutor4)
available to be called from Python (in shared lib ``egsfortran``).

Global variables were defined in Python to match the COMMON blocks
in Fortran, to allow both Python and Fortran to share the same variables.
These definitions were created in separate repo ``egsnrc`` and copied
here.

Initially, ``tutor4.mortran`` and the ``shower`` main loop were translated
to Python code.  Going forward, more pieces of the original Fortran are
being replaced.

So far this is working only for linux; f2py seems to be more difficult
to get working properly under Windows.


Pre-reqs
--------

* Python 3.9 available and callable with ``python3.9``
* EGSnrc installed and config set up for ``linux``
* EGS_HOME set to ``.../egsnrc/egsnrc/egs_home``


Steps
-----

* in ``egsnrc/egsnrc/egs_home/tutor4`` folder:

* ``python3.9 make_egsfortran.py``.  This compiles the mortran code, runs numpy.f2py
  to create the ``egsfortran...so`` file, and copies it to the
  main ``egsnrc`` folder so ``egsfortran`` is importable as part of the package.
* ``python3.9 tutor4.py`` to run the Python version

In the above, there may be isses which are best seen using the standard ``make``.
The tutor4 mortran code here is no longer runnable through the usual EGS run,
if you wish to see the original output you will need to temporarily set ``EGS_HOME``
back to the user_code directory in the EGSnrc installed location, then in
that tutor4 folder, run:

* ``make``
* ``tutor4 -p tutor_data``

``tutor4`` outputs details of the particles and interactions at each step.  This
provides a validation baseline as new Python code is introduced.


Python / egsfortran Tips
------------------------

All variables have been converted to lower case in Python.  Fortran ignores
case in variable names, so care must be taken not to leave any in 
uppercase or mixed-case, which Python will treat as a distinct variable,
leading to subtle bugs.

Fortran has 1-based array indexing by default while
Python is zero-based, so care must be taken to ensure indices 
are modified (usually subtract 1) in Python.  Note, however, that
some EGSnrc arrays override the default and specify a zero-based array.'
In Python code, the difference of 1 in array indexing is handled only when
accessing array elements.  Often a separate variable is created to do
the subtraction only once, and to make it easier to find later if moving to
a Python-only implementation.  E.g. ``medium`` is often used as an
array index.  In Python code, we define ``medium_m1 = medium - 1`` and 
use that for array indexing.

All globals used in the program must be defined also in Python.  If not, 
then setting or reading values will actually be affecting a local
Python variable, and this can be very difficult to debug.  In some cases
this could mean that a Fortran variable has kept its default value (often 0)
rather than what you are trying to set.

Assigning immutables to a global name must use the common block name 
and attribute, e.g. ``stack.np = 1``, otherwise Python will re-bind
the name to a non-egsfortran Python local/global.
This is not specific to this code, but normal Python behavior for
binding of variable names.

Setting elements of arrays can be done without the common block reference,
as this changes the internal state of the array, not where
the variable name is pointing. E.g. ``iq[0] = 1`` is okay without the
common block reference.  *Using* variable values, as part of a right-hand-side
expression, can always be done without concern.
