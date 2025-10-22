import json
import sys

FILE_HEADER = """#pragma once
#include <array>
#include <string_view>

// THIS FILE WAS GENERATED, CHANGES WON'T BE PERSISTED!

// Note that this truncates the DB before inserting!
// We do this to avoid clashes with existing data, but we can
// revisit this later on if needed.
// The MSVC compiler limits the maximum literal string size to 16Kb,
// So we have to split this into chunks

constexpr std::array<std::string_view, 2> create_default_monomers{
    R"SQL(
"""

SQL_INSERT_PREAMBLE = """BEGIN TRANSACTION;
DELETE FROM monomer_definitions;
"""

SQL_INSERT_TEMPLATE = (
    "INSERT INTO monomer_definitions VALUES("
    "{idx},'{symbol}','{polymer_type}','{natural_analog}','{smiles}',"
    "'{core_smiles}','{name}','{monomer_type}','{author}','{pdbcode}');\n")

SQL_COMMIT = "COMMIT;\n"

BLOCK_SEPARATOR = """)SQL",
    R"SQL(
"""

FILE_FOOTER = """)SQL"};
"""


def format_monomer(idx, monomer):
    monomer = {k.lower(): v for k, v in monomer.items()}
    return SQL_INSERT_TEMPLATE.format(idx=idx, **monomer)


def writer(f, chunk_size):
    block_size = 0

    def write(txt):
        nonlocal block_size
        sz = len(txt)
        if block_size + sz >= chunk_size:
            f.write(BLOCK_SEPARATOR)
            block_size = 0
        f.write(txt)
        block_size += sz

    return write


def write_monomer_header(monomers, out_file_name):

    with open(out_file_name, 'wt') as f:
        f.write(FILE_HEADER)

        # MSVC doesn't support literal strings larger than 16kb,
        # so we need to split the SQL into chunks
        write_and_count = writer(f, 16 * 1024)
        write_and_count(SQL_INSERT_PREAMBLE)
        for idx, monomer in enumerate(monomers):
            write_and_count(format_monomer(idx, monomer))
        write_and_count(SQL_COMMIT)

        f.write(FILE_FOOTER)


def main(args):
    if len(args) != 3:
        print(
            f'Usage: {args[0]} [input json definitions] [header output file name]'
        )
        exit(1)

    # Load the monomer definitions
    with open(args[1]) as f:
        monomers = json.load(f)

    write_monomer_header(monomers, args[2])


if __name__ == '__main__':
    main(sys.argv)
