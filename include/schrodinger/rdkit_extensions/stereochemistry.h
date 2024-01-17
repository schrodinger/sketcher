#pragma once

#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"

// Forward declarations:
namespace RDKit
{
class Atom;
class Conformer;
class Bond;
class ROMol;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

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
 * @return the chiral label for the given atom.
 */
RDKIT_EXTENSIONS_API std::string
get_atom_chirality_label(const RDKit::Atom& atom);

/**
 * @return the chiral label for the given bond.
 */
RDKIT_EXTENSIONS_API std::string get_bond_stereo_label(const RDKit::Bond& bond);

/**
 * Set enhanced stereo for any chiral atoms that don't already have it.
 * If the chiral flag is on, ungrouped chiral atoms will go into the ABS
 * group. If the flag is off or not present, ungrouped atoms will go into
 * a new AND group.
 */
RDKIT_EXTENSIONS_API void
add_enhanced_stereo_to_chiral_atoms(RDKit::ROMol& mol);

/**
 * A custom wrapper around RDkit's WedgeMolBonds() that makes sure
 * we don't wedge attachment point dummy atoms
 */
RDKIT_EXTENSIONS_API void wedgeMolBonds(RDKit::ROMol& mol,
                                        const RDKit::Conformer* conf);

} // namespace rdkit_extensions
} // namespace schrodinger