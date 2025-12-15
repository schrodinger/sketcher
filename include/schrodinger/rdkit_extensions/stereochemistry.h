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
 * A custom wrapper around RDkit's WedgeMolBonds() that makes sure
 * we don't wedge attachment point dummy atoms
 */
RDKIT_EXTENSIONS_API void wedgeMolBonds(RDKit::ROMol& mol,
                                        const RDKit::Conformer* conf);

} // namespace rdkit_extensions
} // namespace schrodinger