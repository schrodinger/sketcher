#!/usr/bin/env python3

import re
import subprocess
import sys


def main():
    files = sys.argv[1:]
    if not files:
        return 0

    # Filter for C/C++ files only
    cpp_pattern = re.compile(r'\.(cpp|h)$')
    cpp_files = [f for f in files if cpp_pattern.search(f)]

    if not cpp_files:
        return 0

    bad_files = []
    for filename in cpp_files:
        # Get formatted content from clang-format
        try:
            formatted_content = subprocess.check_output(
                ["clang-format", filename],
                universal_newlines=True
            )
        except subprocess.CalledProcessError:
            print(f"ERROR: clang-format failed on {filename}")
            return 1

        # Read actual file content
        with open(filename, 'r') as fh:
            full_content = fh.read()

        # Compare
        if formatted_content != full_content:
            bad_files.append(filename)

    if bad_files:
        print("ERROR: The following files are non-conforming to clang-format style:")
        print()
        print("\n".join(bad_files))
        print()
        print("To fix formatting, run:")
        print("  clang-format -i " + " ".join(bad_files))
        print()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
