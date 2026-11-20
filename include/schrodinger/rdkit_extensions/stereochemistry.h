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
 * @return whether the given bond is drawn as a wedge or a dash, i.e. whether it
 * has a bond direction of BEGINWEDGE or BEGINDASH
 */
RDKIT_EXTENSIONS_API bool is_wedge_or_dash(const RDKit::Bond& bond);

/**
 * A custom wrapper around RDkit's WedgeMolBonds() that makes sure
 * we don't wedge attachment point dummy atoms
 *
 * @param mol The molecule to wedge
 * @param conf The conformer to use when determining wedge directions
 * @param keep_existing_wedges If false, all existing single bond directions
 * will be cleared and wedges will be recalculated for all chiral centers. If
 * true, existing wedged, dashed, and wavy bonds will be kept, and wedges will
 * only be added for chiral centers that don't already have one. There are two
 * exceptions where an existing wedge or dash won't be kept as-is:
 *   - If it starts at a chiral center and disagrees with that center's chiral
 *     tag, then it's flipped (wedge to dash or vice versa) so that it agrees.
 *     In other words, the chiral tags take precedence over the bond directions,
 *     so stale bond directions can't change the stereochemistry.
 *   - If it starts at a non-chiral atom and ends at a chiral center that has no
 *     wedge or dash of its own, then it's cleared.  Otherwise, RDKit would
 *     consider that chiral center to be wedged already and wouldn't give it a
 *     wedge that describes its stereochemistry.
 * Regardless of this value, any double bond direction markers on single bonds
 * (i.e. "/" and "\" from SMILES) are cleared, as is the MDL wedging info for
 * wedged and dashed bonds. The MDL wedging info for wavy single bonds and
 * crossed double bonds is only cleared if this value is false.
 */
RDKIT_EXTENSIONS_API void wedgeMolBonds(RDKit::ROMol& mol,
                                        const RDKit::Conformer* conf,
                                        bool keep_existing_wedges = false);

} // namespace rdkit_extensions
} // namespace schrodinger