/* -------------------------------------------------------------------------
 * Implements schrodinger::rdkit_extensions:: miscellaneous mol operations
 mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#pragma once

#include <string>
#include <unordered_map>

#include "schrodinger/rdkit_extensions/molops.h"

#include <rdkit/GraphMol/MolOps.h>

namespace schrodinger
{
namespace rdkit_extensions
{

const std::string ATTACHMENT_POINT_LABEL_PREFIX{"_AP"};
const std::string R_GROUP_LABEL_PREFIX{"_R"};
const int DUMMY_ATOMIC_NUMBER{0};

// Map from the one-letter code of each of the 20 standard amino acids to its
// upper-case three-letter (PDB residue name) code
const std::unordered_map<std::string, std::string>
    AMINO_ACID_THREE_LETTER_CODES({
        {"A", "ALA"}, // Alanine
        {"R", "ARG"}, // Arginine
        {"N", "ASN"}, // Asparagine
        {"D", "ASP"}, // Aspartic
        {"C", "CYS"}, // Cysteine
        {"Q", "GLN"}, // Glutamine
        {"E", "GLU"}, // Glutamic
        {"G", "GLY"}, // Glycine
        {"H", "HIS"}, // Histidine
        {"I", "ILE"}, // Isoleucine
        {"L", "LEU"}, // Leucine
        {"K", "LYS"}, // Lysine
        {"M", "MET"}, // Methionine
        {"F", "PHE"}, // Phenylalanine
        {"P", "PRO"}, // Proline
        {"S", "SER"}, // Serine
        {"T", "THR"}, // Threonine
        {"W", "TRP"}, // Tryptophan
        {"Y", "TYR"}, // Tyrosine
        {"V", "VAL"}, // Valine
    });

} // namespace rdkit_extensions
} // namespace schrodinger
