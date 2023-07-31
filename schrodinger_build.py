"""
Convenience script for building/testing

> source $SCHRODINGER_SRC/mmshare/build_env
> $SCHRODINGER/run schrodinger_build.py
"""

import argparse
import os
import subprocess
import sys

import pytest
from schrodinger.test import add_build_tools_to_path

with add_build_tools_to_path():
    from library_definitions import get_library_definitions
    LIBDEFS = get_library_definitions()

PREFIX_PATHS = (
    LIBDEFS["QT"].base_directory,
    # LIBDEFS["RDKIT"].base_directory,
    # FIXME: lib/cmake/rdkit/rdkit-targets.cmake required updating these paths:
    # RDKit::rdkit_base PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
    # RDKit::rdkit_py_base PROPERTIES INTERFACE_LINK_LIBRARIES
    # RDKit::Depictor PROPERTIES INTERFACE_LINK_LIBRARIES
    # RDKit::FileParsers PROPERTIES INTERFACE_LINK_LIBRARIES
    # RDKit::MolDraw2D PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
    # RDKit::MolDraw2D PROPERTIES INTERFACE_LINK_LIBRARIES
    "/Users/vonbarge/schrodinger/sketcher_libcopy/rdkit-Release_2023_03_2-BLDMGR-7314",
    # FIXME: 3.4.0 doesn't have cmake files!
    LIBDEFS["EIGEN"].base_directory + "/../eigen-0eeb60526ce0/share/eigen3",
    LIBDEFS["COORDGENLIBS"].library_directory,  # for maeparser
    LIBDEFS["BOOST"].base_directory,
    LIBDEFS["FMTLIB"].base_directory,
)

CONFIG_FLAGS = [
    f"-DCMAKE_PREFIX_PATH={';'.join(PREFIX_PATHS)}",
    "-DCMAKE_OSX_ARCHITECTURES=x86_64",  # SCHRODINGER_LIB uses x86_64 on Darwin
]


def parse_args(argv=None):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    msg = "Do a build from scratch."
    parser.add_argument("--clean", action="store_true", help=msg)

    msg = "Run tests after building."
    parser.add_argument("--with-tests", action="store_true", help=msg)

    msg = "Pass the given arguments to pytest."
    parser.add_argument("--test-args", dest="TEST_ARGS", default="", help=msg)

    return parser.parse_args(argv)


def main(argv=None):

    args = parse_args(argv)
    top_dir = os.path.dirname(os.path.abspath(sys.argv[0]))
    build_dir = os.path.join(top_dir, "build")

    run_cmake = lambda x: subprocess.check_call(["cmake"] + x)

    if args.clean:
        run_cmake(["--build", build_dir, "--target", "clean"])
    run_cmake([top_dir, "-B", build_dir] + CONFIG_FLAGS)
    run_cmake(["--build", build_dir])

    if args.with_tests:
        pytest_args = [str(build_dir)] + args.TEST_ARGS.split()
        pytest.main(pytest_args)


if __name__ == "__main__":
    main()
