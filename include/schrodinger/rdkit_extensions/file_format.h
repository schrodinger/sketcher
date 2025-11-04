/* -------------------------------------------------------------------------
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <sstream>
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
    // Atomistic-specific formats
    SMILES,
    EXTENDED_SMILES,
    SMARTS,
    EXTENDED_SMARTS,
    MDL_MOLV2000,
    MDL_MOLV3000,
    MAESTRO,
    INCHI,
    INCHI_KEY, // only used to write
    PDB,
    MOL2,
    XYZ,
    MRV,
    CDXML,
    // Monomer-specific formats
    HELM,
    FASTA_PEPTIDE, // only used to read
    FASTA_DNA,     // only used to read
    FASTA_RNA,     // only used to read
    FASTA,         // only used to write
};

/**
 * Enum class specifying supported compression types
 *
 * 'UNKNOWN' means the file is either not compressed, or uses an unknown
 * compression method.
 */
enum class RDKIT_EXTENSIONS_API CompressionType { UNKNOWN, GZIP, ZSTD };

RDKIT_EXTENSIONS_API std::ostream&
operator<<(std::ostream& os, const CompressionType& compression_type);

/**
 * All supported formats for standard atomistic molecules, monomeric
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
RDKIT_EXTENSIONS_API std::vector<std::string>
get_mol_extensions(const Format format);
RDKIT_EXTENSIONS_API std::vector<std::string>
get_rxn_extensions(const Format format);
RDKIT_EXTENSIONS_API std::vector<std::string>
get_seq_extensions(const Format format);

/**
 * @param filename file path
 * @return corresponding file format enum
 * @throw std::invalid_argument if file extension is unrecognized
 */
RDKIT_EXTENSIONS_API Format
get_file_format(const boost::filesystem::path& filename);

/**
 * Determines the type of the compression applied to a given file
 * by reading the magic number from the beginning of the file.
 *
 * @param filename file path
 * @return corresponding compression type enum
 * @throw std::invalid_argument if file cannot be opened
 */
[[nodiscard]] RDKIT_EXTENSIONS_API CompressionType
get_compression_type(const boost::filesystem::path& filename);

/**
 * Determines the type of the compression applied to a string buffer
 * by reading the magic number from the beginning of the string
 *
 * @param a string stream
 * @return corresponding compression type enum
 * @throw std::invalid_argument if string stream is invalid
 */
[[nodiscard]] RDKIT_EXTENSIONS_API CompressionType
get_compression_type(std::istringstream& sstream);

/**
 * Determines the type of the compression applied to a given file
 * by looking at the extension
 *
 * @param filename file path
 * @return corresponding compression type enum
 */
[[nodiscard]] RDKIT_EXTENSIONS_API CompressionType
get_compression_type_from_ext(const boost::filesystem::path& filename);

} // namespace rdkit_extensions
} // namespace schrodinger
