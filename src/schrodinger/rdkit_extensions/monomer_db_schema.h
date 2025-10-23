#pragma once
#include <string_view>

constexpr std::string_view create_monomer_definitions_table_sql = R"SQL(
BEGIN TRANSACTION;
CREATE TABLE IF NOT EXISTS "monomer_definitions" (
    "ID"              INTEGER,
    "SYMBOL"          TEXT,
    "POLYMER_TYPE"    TEXT,
    "NATURAL_ANALOG"  TEXT,
    "SMILES"          TEXT UNIQUE,
    "CORE_SMILES"     TEXT,
    "NAME"            TEXT UNIQUE,
    "MONOMER_TYPE"    TEXT,
    "AUTHOR"          TEXT,
    "PDBCODE"         TEXT,
    PRIMARY KEY("ID" AUTOINCREMENT)
);
CREATE UNIQUE INDEX "monomertype_symbol_idx" ON "monomer_definitions" (
    "POLYMER_TYPE"    ASC,
    "SYMBOL"          ASC
);
CREATE INDEX "pdbcode_idx" on "monomer_definitions" ( "PDBCODE" ASC );
CREATE INDEX "core_smiles_idx" on "monomer_definitions" ( "CORE_SMILES" ASC );
COMMIT;
)SQL";