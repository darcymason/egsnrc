#!/usr/bin/env python

import sys
from setuptools import setup, find_packages
from pathlib import Path


description = "egsnrc is the EGSnrc Monte Carlo code ported to Python"

needs_pytest = {"pytest", "test", "ptr"}.intersection(sys.argv)
pytest_runner = ["pytest-runner"] if needs_pytest else []

_py_modules = []

CLASSIFIERS = [
    "License :: OSI Approved :: GNU General Public License (GPL)",
    "Intended Audience :: Healthcare Industry",
    "Intended Audience :: Science/Research",
    # "Development Status :: 5 - Production/Stable",
    "Programming Language :: Python :: 3.8",
    "Programming Language :: Python :: 3.9",
    "Operating System :: OS Independent",
    "Topic :: Scientific/Engineering :: Medical Science Apps.",
    "Topic :: Scientific/Engineering :: Physics",
]

KEYWORDS = "python ionizing radiation simulation transpiler"

URL = "https://github.com/darcymason/egsnrc2py/egsnrc"
REQUIRES = []
SETUP_REQUIRES = pytest_runner

# get long description from README.md
BASE_PATH = Path(__file__).resolve().parent
with open(BASE_PATH / "README.md", 'r') as f:
    LONG_DESCRIPTION = f.read()

# Get __version from _version.py
with open(BASE_PATH / "egsnrc" / "_version.py", 'r') as f:
    exec(f.read())


opts = dict(
    name="egsnrc",
    python_requires='>=3.8',
    version=__version__,  # noqa
    maintainer="Darcy Mason and contributors",
    maintainer_email="darcymason@gmail.com",
    author="Darcy Mason and contributors",
    author_email="darcymason@gmail.com",
    description=description,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=URL,
    # download_url=DOWNLOAD_URL,
    # license=LICENSE"GPLv3",
    keywords=KEYWORDS,
    classifiers=CLASSIFIERS,
    # packages=find_packages(),
    py_modules=_py_modules,
    # package_data=PACKAGE_DATA,
    # include_package_data=True,
    install_requires=REQUIRES,
    setup_requires=SETUP_REQUIRES,
    tests_require=['pytest'],
    zip_safe=False,
)

if __name__ == "__main__":
    setup(**opts)
