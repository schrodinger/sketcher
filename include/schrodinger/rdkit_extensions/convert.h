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
class Atom;
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
    MDL_MOLV2000,
    MDL_MOLV3000,
    INCHI,
    INCHI_KEY,
    PDB,
};

/**
 * @param text input text block
 * @param format specified format from which to interpret the text
 * @return resulting rdkit molecule
 * @throw std::invalid_argument if text is not valid for the given format
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
text_to_rdmol(const std::string& text, Format format = Format::AUTO_DETECT);

/**
 * @param text input text block
 * @param format specified format from which to interpret the text
 * @return resulting rdkit reaction
 * @throw std::invalid_argument if text is not valid for the given format
 */
RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::ChemicalReaction>
text_to_reaction(const std::string& text, Format format = Format::AUTO_DETECT);

/**
 * @param mol rdkit molecule
 * @param format specified format to which to serialize
 * @return requested text representation of the given molecule
 */
RDKIT_EXTENSIONS_API std::string rdmol_to_text(const RDKit::ROMol& mol,
                                               Format format);

/**
 * @param mol rdkit molecule
 * @param format specified format to which to serialize
 * @return requested text representation of the given reaction
 */
RDKIT_EXTENSIONS_API std::string
reaction_to_text(const RDKit::ChemicalReaction& reaction, Format format);

/**
 * @param atom rdkit atom
 * @return whether the atom is an attachment point
 */
RDKIT_EXTENSIONS_API bool is_attachment_point_dummy(const RDKit::Atom& atom);

/**
 * Set enhanced stereo for any chiral atoms that don't already have it.
 * If the chiral flag is on, ungrouped chiral atoms will go into the ABS
 * group. If the flag is off or not present, ungrouped atoms will go into
 * a new AND group.
 */
RDKIT_EXTENSIONS_API void
add_enhanced_stereo_to_chiral_atoms(RDKit::ROMol& mol);

/**
 * Set bond directions based on the input data from SDF
 * @param rdk_mol rdkit mol
 */
RDKIT_EXTENSIONS_API void reapply_molblock_wedging(RDKit::ROMol& rdk_mol);

} // namespace rdkit_extensions
} // namespace schrodinger
