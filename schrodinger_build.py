"""
Convenience script for building/testing in the Schrodinger environment

> source $SCHRODINGER_SRC/mmshare/build_env
> $SCHRODINGER/run schrodinger_build.py
"""

import argparse
import os
import subprocess
import sys

from schrodinger.test import add_build_tools_to_path

with add_build_tools_to_path():
    from library_definitions import get_library_definitions
    LIBDEFS = get_library_definitions()

PREFIX_PATHS = (
    LIBDEFS["BOOST"].base_directory,
    LIBDEFS["FMTLIB"].base_directory,
    LIBDEFS["ZSTD"].base_directory,
    LIBDEFS["RDKIT"].base_directory,
    LIBDEFS["QT"].base_directory,
)

CONFIG_FLAGS = [
    f"-DCMAKE_PREFIX_PATH={';'.join(PREFIX_PATHS)}",
    "-DCMAKE_OSX_ARCHITECTURES=x86_64",  # SCHRODINGER_LIB uses x86_64
    "-DBUILD_SHARED_LIBS=OFF",
    "-DCMAKE_BUILD_TYPE=Release",
]


def parse_args(argv=None):

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    msg = "Run a clean build"
    parser.add_argument("--clean", action="store_true", help=msg)

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


if __name__ == "__main__":
    main()
