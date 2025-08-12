#include "schrodinger/sketcher/molviewer/monomer_utils.h"

#include <algorithm>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/AtomIterators.h>
#include <rdkit/GraphMol/ROMol.h>

namespace
{
const std::string MONOMERIC_PROPERTY = "SKETCHER_ATOM_IS_MONOMERIC";
}

namespace schrodinger
{
namespace sketcher
{

void set_atom_monomeric(RDKit::Atom* const atom)
{
    atom->setProp(MONOMERIC_PROPERTY, true);
}

void clear_monomeric_property(RDKit::Atom* const atom)
{
    atom->clearProp(MONOMERIC_PROPERTY);
}

bool is_atom_monomeric(const RDKit::Atom* const atom)
{
    bool is_monomeric = false;
    atom->getPropIfPresent(MONOMERIC_PROPERTY, is_monomeric);
    return is_monomeric;
}

bool contains_monomeric_atom(const RDKit::ROMol& mol)
{
    auto atoms_iter = mol.atoms();
    return std::any_of(atoms_iter.begin(), atoms_iter.end(), is_atom_monomeric);
}

} // namespace sketcher
} // namespace schrodinger
