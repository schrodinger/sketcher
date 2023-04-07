/* -------------------------------------------------------------------------
 * Declares schrodinger::rdkit_extensions:: text block <-> rdkit mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <string>

#include <boost/shared_ptr.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"

// Forward declarations:
namespace RDKit
{
class ROMol;
class RWMol;
class ChemicalReaction;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

/**
 * Available input/output formats for import/export
 */
enum class Format {
    AUTO_DETECT,
    SMILES,
    EXTENDED_SMILES,
    SMARTS,
    MAESTRO,
    MDL_MOLV2000,
    MDL_MOLV3000,
    INCHI,
    INCHI_KEY,
    PDB,
    RDMOL_BINARY_BASE64,
    XYZ,
};

/**
 * @param text input text block
 * @param format specified format from which to interpret the text
 * @return resulting rdkit molecule
 * @throw std::invalid_argument if text is not valid for the given format
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
to_rdkit(const std::string& text, Format format = Format::AUTO_DETECT);

/**
 * @param text input text block
 * @param format specified format from which to interpret the text
 * @return resulting rdkit reaction
 * @throw std::invalid_argument if text is not valid for the given format
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::ChemicalReaction>
to_rdkit_reaction(const std::string& text, Format format = Format::AUTO_DETECT);

/**
 * @param mol rdkit molecule
 * @param format specified format to which to serialize
 * @return requested text representation of the given molecule
 */
RDKIT_EXTENSIONS_API std::string to_string(const RDKit::ROMol& mol,
                                           Format format);

/**
 * @param mol rdkit molecule
 * @param format specified format to which to serialize
 * @return requested text representation of the given reaction
 */
RDKIT_EXTENSIONS_API std::string
to_string(const RDKit::ChemicalReaction& reaction, Format format);

} // namespace rdkit_extensions
} // namespace schrodinger
