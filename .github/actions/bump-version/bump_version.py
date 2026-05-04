#!/usr/bin/env python3
"""Compute the next sketcher version from the current one.

Default mode increments the build number (third field).
With ``--bump-quarter``, rolls year/quarter forward and resets build to 0
(Q4 wraps to Q1 of the next year).
"""
import argparse
import sys


def next_version(current: str, bump_quarter: bool = False) -> str:
    """Return the version that should follow ``current``.

    ``current`` is expected in ``YYYY.Q.B`` form.
    Raises ``ValueError`` on malformed input or out-of-range quarter.
    """
    parts = current.split(".")
    if len(parts) != 3:
        raise ValueError(
            f"expected version in YYYY.Q.B form, got {current!r}")
    try:
        year, quarter, build = (int(p) for p in parts)
    except ValueError as err:
        raise ValueError(
            f"version components must be integers, got {current!r}") from err
    if not 1 <= quarter <= 4:
        raise ValueError(
            f"quarter must be in 1..4, got {quarter} (from {current!r})")
    if build < 0:
        raise ValueError(
            f"build must be non-negative, got {build} (from {current!r})")
    if bump_quarter:
        if quarter == 4:
            return f"{year + 1}.1.0"
        return f"{year}.{quarter + 1}.0"
    return f"{year}.{quarter}.{build + 1}"


def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute the next sketcher version.")
    parser.add_argument("version", help="current version string (YYYY.Q.B)")
    parser.add_argument(
        "--bump-quarter",
        action="store_true",
        help="roll year/quarter forward and reset build to 0")
    args = parser.parse_args(argv)
    try:
        print(next_version(args.version, bump_quarter=args.bump_quarter))
    except ValueError as err:
        print(f"error: {err}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
