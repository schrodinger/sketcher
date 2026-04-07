#include "schrodinger/rdkit_extensions/monomer_utils.h"

#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include <memory>

namespace schrodinger
{
namespace rdkit_extensions
{

/// Return a map from atom index to attachment point number.
std::vector<unsigned int> make_attch_map(const RDKit::ROMol& mol)
{
    std::vector<unsigned int> attch_map(mol.getNumAtoms(), NO_ATTACHMENT);
    for (const auto atom : mol.atoms()) {
        if (atom->hasProp(RDKit::common_properties::molAtomMapNumber)) {
            attch_map[atom->getIdx()] = atom->getProp<unsigned int>(
                RDKit::common_properties::molAtomMapNumber);
        }
    }
    return attch_map;
}

void neutralizeAtoms(RDKit::ROMol& mol)
{
    // Algorithm for neutralizing molecules from
    // https://www.rdkit.org/docs/Cookbook.html#neutralizing-molecules by Noel
    // O’Boyle Will neutralize the molecule by adding or removing hydrogens as
    // needed. This will ensure SMILES can be used to match atomistic structures
    // to the correct monomer.
    static const std::unique_ptr<RDKit::RWMol> neutralize_query(
        RDKit::SmartsToMol(
            "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"));
    for (const auto& match : RDKit::SubstructMatch(mol, *neutralize_query)) {
        auto atom = mol.getAtomWithIdx(match[0].second);
        auto chg = atom->getFormalCharge();
        auto hcount = atom->getTotalNumHs();
        atom->setFormalCharge(0);
        atom->setNumExplicitHs(hcount - chg);
        atom->updatePropertyCache();
    }

    // This query matches the dipole representation of sulfoxides, nitro, etc.
    static const std::unique_ptr<RDKit::RWMol> dipole_query(
        RDKit::SmartsToMol("[S+]-[O-]"));
    for (const auto& match : RDKit::SubstructMatch(mol, *dipole_query)) {
        auto a1 = mol.getAtomWithIdx(match[0].second);
        auto a2 = mol.getAtomWithIdx(match[1].second);
        a1->setFormalCharge(0);
        a2->setFormalCharge(0);
        auto b = mol.getBondBetweenAtoms(a1->getIdx(), a2->getIdx());
        b->setBondType(RDKit::Bond::DOUBLE);
        a1->updatePropertyCache();
        a2->updatePropertyCache();
        b->updatePropertyCache();
    }
}

} // namespace rdkit_extensions
} // namespace schrodinger
