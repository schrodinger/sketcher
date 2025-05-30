"""
A unit test to ensure that the versions of dependency packages stay in sync
between open-source builds of Sketcher and the builds used in the Schrodinger
core suite.

Note that this test should **not** be launched from premake or the generated
Makefile. Buildinger strips $SCHRODINGER_PACKAGE_FACTORY_BASE from the
environment when it executes those, which means that the unit test may not be
able to find the pixi executable, which it needs to fetch version information.
The test should instead be launched directly from a terminal via::

    $SCHRODINGER/utilities/py.test test_dependency_versions.py
"""

import itertools
import json
import logging
import os
import pathlib
import re
import sys
import types
from collections.abc import Iterable

SKETCHER_TEST_DIR = pathlib.Path(__file__).parent
SKETCHER_VERSION_FILE = SKETCHER_TEST_DIR / ".." / "external" / "versions.json"
# The name of the package factory virtual environment to use when fetching
# version info from mmshare
MMSHARE_PF_ENVIRONMENT = "schrodinger"
# Packages that we should exclude from the version checks
PACKAGES_TO_IGNORE = [
    # Core suite builds frequently use a commercial-only LTS release of Qt,
    # which are not available for open-source builds of Sketcher
    "qt",
    # The core suite doesn't include emscripten
    "emscripten",
    # this repo is currently on a newer RDKit and Boost versions than the core
    # suite
    "rdkit",
    "boost",
    # fmt 10.1.1 reports itself as version 10.1.0, which confuses CMake, so we
    # use 10.1.0 for this repo
    "fmt"
]


def import_package_factory() -> types.ModuleType:
    """
    Import the package_factory module from the mmshare repo
    """
    schrodinger_src = pathlib.Path(os.environ["SCHRODINGER_SRC"])
    build_tools_path = schrodinger_src / 'mmshare' / 'build_tools'
    if not build_tools_path.is_dir():
        raise RuntimeError(
            "$SCHRODINGER_SRC/mmshare/build_tools cannot be found.")
    sys.path.insert(0, str(build_tools_path))
    try:
        import package_factory
    except Exception as err:
        raise RuntimeError("The package factory module could not be imported "
                           f"from {build_tools_path}:\n"
                           f"\t{str(err)}")
    finally:
        sys.path.pop(0)
    return package_factory


package_factory = import_package_factory()


def get_dependency_versions_from_mmshare(
        packages: Iterable) -> dict[str, tuple[int]]:
    """
    Get information about package versions from mmshare.

    :return: Version information as a dictionary of
      {package name: (major version, minor version, patch version)}
    """
    logger = logging.getLogger("test_dependency_versions")
    mmshare_pf = package_factory.Pixi(logger, MMSHARE_PF_ENVIRONMENT)
    mmshare_manifest_file = package_factory.Pixi.get_requirements_file()
    # Build a pattern that matches only the packages we're looking for. The
    # carat and dollar sign specify that the package name must be an exact match
    # without any extra preceding or trailing characters, since we don't want to
    # match, e.g., xgboost when looking for boost.
    packages_description = "|".join(f"^{dep}$" for dep in packages)
    # --no-lockfile-update tells Pixi not to make any changes to the mmshare
    # build env, even if it's out of date or not built yet
    cmd = [
        "list", packages_description, "--manifest-path",
        str(mmshare_manifest_file), "-e", MMSHARE_PF_ENVIRONMENT,
        "--no-lockfile-update", "--json"
    ]
    json_data = mmshare_pf.run(cmd, capture_output=True)
    dep_info = json.loads(json_data)
    return {
        cur_dep["name"]: normalize_version(cur_dep['version'])
        for cur_dep in dep_info
    }


def normalize_version(version_text: str) -> tuple[int]:
    """
    Normalize the version strings that have been retrieved from package factory
    and convert everything to tuples of ints. For example, this method will
    convert
      - 1_81_0_py311_cpp20
      - 3.4.0_e67c494c_bldmgr_9082
      - v1.5.5_boost_1_81_0_BLDMGR_7247
      - release_2024_9_6_bldmgr_9464
    to (1, 81, 0), (3, 4, 0), (1, 5, 5), and (2024, 9, 6)
    """
    # remove any leading "v" or "release" prefixes
    version_text = version_text.removeprefix("v")
    version_text = version_text.removeprefix("release")
    version_text = version_text.lstrip("_")
    # ignore any fields that contains letters, as well as anything afterwards
    version_fields = re.split(r"[._]", version_text)
    numeric_fields = itertools.takewhile(lambda field: field.isdigit(),
                                         version_fields)
    return tuple(map(int, numeric_fields))


def test_dependency_versions():
    """
    Ensure that all dependency versions are in sync (other than packages listed
    in PACKAGES_TO_IGNORE).
    """
    with open(SKETCHER_VERSION_FILE) as sketcher_versions_in:
        sketcher_versions = json.load(sketcher_versions_in)
    for package in PACKAGES_TO_IGNORE:
        sketcher_versions.pop(package)
    # convert the versions to tuples of ints.  Otherwise, RDKit version
    # comparisons may fail because of leading zeroes in the month.  (I.e.
    # "2024.09.6" != "2024.9.6", even though they represent the same versions)
    sketcher_versions_as_ints = {
        package: tuple(map(int, version.split(".")))
        for package, version in sketcher_versions.items()
    }
    mmshare_versions = get_dependency_versions_from_mmshare(
        sketcher_versions.keys())
    assert sketcher_versions_as_ints == mmshare_versions
