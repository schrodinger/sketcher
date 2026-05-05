"""Unit tests for the bump_version helper used by the bump-version GHA."""
import pathlib
import sys

import pytest

SKETCHER_REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
BUMP_VERSION_DIR = (SKETCHER_REPO_ROOT / ".github" / "actions" / "bump-version")
sys.path.insert(0, str(BUMP_VERSION_DIR))
try:
    import bump_version  # noqa: E402
finally:
    sys.path.pop(0)


@pytest.mark.parametrize("current,expected", [
    ("2026.3.36", "2026.3.37"),
    ("2026.3.0", "2026.3.1"),
    ("2027.1.99", "2027.1.100"),
    ("2026.4.7", "2026.4.8"),
])
def test_default_increments_build(current, expected):
    assert bump_version.next_version(current) == expected


@pytest.mark.parametrize("current,expected", [
    ("2026.1.5", "2026.2.0"),
    ("2026.2.5", "2026.3.0"),
    ("2026.3.36", "2026.4.0"),
    ("2026.4.5", "2027.1.0"),
    ("2026.4.0", "2027.1.0"),
    ("2099.4.999", "2100.1.0"),
])
def test_bump_quarter(current, expected):
    assert bump_version.next_version(current, bump_quarter=True) == expected


@pytest.mark.parametrize("malformed", [
    "",
    "2026.3",
    "2026.3.36.1",
    "abc.3.36",
    "2026.x.36",
    "2026.3.x",
    "2026..36",
])
def test_malformed_input_raises(malformed):
    with pytest.raises(ValueError):
        bump_version.next_version(malformed)


@pytest.mark.parametrize("bad_quarter", [
    "2026.0.5",
    "2026.5.5",
    "2026.-1.5",
])
def test_out_of_range_quarter_raises(bad_quarter):
    with pytest.raises(ValueError):
        bump_version.next_version(bad_quarter)


def test_negative_build_raises():
    with pytest.raises(ValueError):
        bump_version.next_version("2026.3.-1")
