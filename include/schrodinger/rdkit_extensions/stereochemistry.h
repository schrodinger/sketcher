#pragma once

#include <optional>
#include <string>

#include <boost/noncopyable.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"

// Forward declarations:
namespace RDKit
{
class Atom;
class Conformer;
class Bond;
class ROMol;
enum class StereoGroupType;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

const std::string HAIR_SPACE = " ";
const std::string ABSOLUTE_STEREO_PREFIX = "abs" + HAIR_SPACE;
const std::string OR_STEREO_PREFIX = "or" + HAIR_SPACE;
const std::string AND_STEREO_PREFIX = "and" + HAIR_SPACE;

using EnhancedStereo = std::pair<RDKit::StereoGroupType, unsigned int>;

class RDKIT_EXTENSIONS_API UseModernStereoPerception : public boost::noncopyable
{
  public:
    UseModernStereoPerception();
    ~UseModernStereoPerception();

  private:
    bool m_stereo_algo_state;
};

/*
 * Assigns stereochemistry to the given molecule, from 3D or 2D if possible
 * @param rdk_mol rdkit mol
 */
RDKIT_EXTENSIONS_API void assign_stereochemistry(RDKit::ROMol& mol);

/**
 * @param strip_abs if true, the absolute stereo label is stripped from the
 * "abs" label
 * @return the chiral label for the given atom.
 */
RDKIT_EXTENSIONS_API std::string
get_atom_chirality_label(const RDKit::Atom& atom, bool strip_abs = false);

/**
 * @return the simplified stereo annotation for the given molecule if the
 * molecule has only one enhanced stereo group that is either an 'or' or 'and'
 * stereo group, otherwise return an empty string.
 */
RDKIT_EXTENSIONS_API std::string
get_simplified_stereo_annotation(const RDKit::ROMol& mol);

/**
 * @return the chiral label for the given bond.
 */
RDKIT_EXTENSIONS_API std::string get_bond_stereo_label(const RDKit::Bond& bond);

/**
 * A custom wrapper around RDkit's WedgeMolBonds() that makes sure
 * we don't wedge attachment point dummy atoms
 */
RDKIT_EXTENSIONS_API void wedgeMolBonds(RDKit::ROMol& mol,
                                        const RDKit::Conformer* conf);

/**
 * Get the enhanced stereo value for the specified atom
 * @return The enhanced stereo type and group write ID for the given atom. (Note
 * that the group ID is only valid for AND or OR stereochemistry, as absolute
 * stereochemistry does not have a group id.) If the atom is not part of any
 * enhanced stereo group, then std::nullopt will be returned.
 */
RDKIT_EXTENSIONS_API std::optional<EnhancedStereo>
get_enhanced_stereo_for_atom(const RDKit::Atom* atom);

/**
 * Modify the enhanced stereo group for the specified atom
 * @param atom The atom to modify
 * @param enh_stereo A pair of the enhanced stereo type and group ID. For AND
 * and OR stereochemistry, this group ID will be set for both the read and write
 * ID. For absolute stereochemistry, this group ID will be ignored.
 */
RDKIT_EXTENSIONS_API void
set_enhanced_stereo_for_atom(RDKit::Atom* atom,
                             const EnhancedStereo& enh_stereo);

} // namespace rdkit_extensions
} // namespace schrodinger