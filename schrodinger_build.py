"""
Convenience script for building/testing
"""

import argparse
import os
import subprocess
from pathlib import Path

import pytest
from schrodinger.test import add_build_tools_to_path

with add_build_tools_to_path():
    from library_definitions import get_library_definitions
    LIBDEFS = get_library_definitions()

LIB_DIR_DICT = {
    "QT": "qt_DIR",
    "RDKIT": "rdkit_DIR",
    "COORDGENLIBS": "maeparser_DIR",
    "BOOST": "boost_DIR",
    "FMTLIB": "fmt_DIR",
}


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

    for lib_name, env_var in LIB_DIR_DICT.items():
        os.putenv(env_var, LIBDEFS[lib_name].base_directory)

    args = parse_args(argv)
    top_dir = Path.cwd()
    build_dir = top_dir / "build"

    run_cmake = lambda x: subprocess.check_call(["cmake"] + x)

    if args.clean:
        run_cmake(["--build", build_dir, "--target", "clean"])
    run_cmake([top_dir, "-B", build_dir])
    run_cmake(["--build", build_dir])

    if args.with_tests:
        pytest_args = [str(build_dir)] + args.TEST_ARGS.split()
        pytest.main(pytest_args)


if __name__ == "__main__":
    main()
