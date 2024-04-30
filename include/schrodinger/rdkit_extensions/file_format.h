/* -------------------------------------------------------------------------
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <string>
#include <vector>

#include <boost/filesystem.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * Enum class specifying supported file formats
 */
enum class Format {
    AUTO_DETECT, // only used to read
    RDMOL_BINARY_BASE64,
    SMILES,
    EXTENDED_SMILES,
    SMARTS,
    MDL_MOLV2000,
    MDL_MOLV3000,
    MAESTRO,
    INCHI,
    INCHI_KEY, // only used to write
    PDB,
    XYZ,
    HELM,
    FASTA_PEPTIDE, // only used to read
    FASTA_DNA,     // only used to read
    FASTA_RNA,     // only used to read
    FASTA,         // only used to write
};

/**
 * All supported formats for standard atomistic molecules, monomeristic
 * sequences, and chemical reactions
 */
extern RDKIT_EXTENSIONS_API const std::vector<Format> MOL_FORMATS;
extern RDKIT_EXTENSIONS_API const std::vector<Format> RXN_FORMATS;
extern RDKIT_EXTENSIONS_API const std::vector<Format> SEQ_FORMATS;

/**
 * @param format file format to consider
 * @return supported molecule/reaction/sequence extensions of the given format
 * @throw std::out_of_range if format is not supported
 */
RDKIT_EXTENSIONS_API std::vector<std::string> get_mol_extensions(Format format);
RDKIT_EXTENSIONS_API std::vector<std::string> get_rxn_extensions(Format format);
RDKIT_EXTENSIONS_API std::vector<std::string> get_seq_extensions(Format format);

/**
 * @param filename file path
 * @return corresponding file format enum
 * @throw std::invalid_argument if file extension is unrecognized
 */
RDKIT_EXTENSIONS_API Format
get_file_format(const boost::filesystem::path& filename);

} // namespace rdkit_extensions
} // namespace schrodinger
