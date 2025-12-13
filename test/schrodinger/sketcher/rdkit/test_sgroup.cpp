
#define BOOST_TEST_MODULE test_sgroup

#include <string>
#include <unordered_set>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/rdkit/sgroup.h"

namespace schrodinger
{
namespace sketcher
{

const std::string S_GROUP_SMILES = "CCCCCCCCC |Sg:n:3,4,5,6:1:ht:::|";
const std::string MULTIPLE_S_GROUP_SMILES =
    "CCCCCCCCCCC "
    "|,,,Sg:n:2,3,4:1:ht:::,Sg:n:7:2:ht:::,Sg:n:1,2,3,4,5,6,7,8:3:ht:::|";

/**
 * Make sure that get_existing_sgroup_for_atoms returns the expected S-groups
 */
BOOST_AUTO_TEST_CASE(test_get_existing_sgroup_for_atoms)
{
    auto mol = rdkit_extensions::to_rdkit(MULTIPLE_S_GROUP_SMILES);
    std::unordered_set<const RDKit::Atom*> atoms({mol->getAtomWithIdx(2),
                                                  mol->getAtomWithIdx(3),
                                                  mol->getAtomWithIdx(4)});
    auto* sgroup = get_existing_sgroup_for_atoms(atoms, *mol);
    BOOST_TEST(sgroup != nullptr);
    BOOST_TEST(get_polymer_label(*sgroup) == "1");

    // remove an atom so that the set no longer matches the S-group
    atoms.erase(atoms.begin());
    sgroup = get_existing_sgroup_for_atoms(atoms, *mol);
    BOOST_TEST(sgroup == nullptr);

    std::unordered_set<const RDKit::Atom*> atoms2({mol->getAtomWithIdx(7)});
    sgroup = get_existing_sgroup_for_atoms(atoms2, *mol);
    BOOST_TEST(sgroup != nullptr);
    BOOST_TEST(get_polymer_label(*sgroup) == "2");

    std::unordered_set<const RDKit::Atom*> atoms3(
        {mol->getAtomWithIdx(1), mol->getAtomWithIdx(2), mol->getAtomWithIdx(3),
         mol->getAtomWithIdx(4), mol->getAtomWithIdx(5), mol->getAtomWithIdx(6),
         mol->getAtomWithIdx(7), mol->getAtomWithIdx(8)});
    sgroup = get_existing_sgroup_for_atoms(atoms3, *mol);
    BOOST_TEST(sgroup != nullptr);
    BOOST_TEST(get_polymer_label(*sgroup) == "3");

    // add another atom so that the set no longer matches the S-group
    atoms3.insert(mol->getAtomWithIdx(9));
    sgroup = get_existing_sgroup_for_atoms(atoms3, *mol);
    BOOST_TEST(sgroup == nullptr);

    // pass in an empty set
    sgroup = get_existing_sgroup_for_atoms({}, *mol);
    BOOST_TEST(sgroup == nullptr);
}

/**
 * Ensure that can_atoms_form_sgroup returns the expected result
 */
BOOST_AUTO_TEST_CASE(test_can_atoms_form_sgroup)
{
    auto mol = rdkit_extensions::to_rdkit("CCCCCCCC");
    std::unordered_set<const RDKit::Atom*> atoms({mol->getAtomWithIdx(2),
                                                  mol->getAtomWithIdx(3),
                                                  mol->getAtomWithIdx(4)});
    // these atoms are contiguous with two connecting bonds, so they can form an
    // S-group
    BOOST_TEST(can_atoms_form_sgroup(atoms, *mol));
    // now the atoms are discontiguous, so they can't
    atoms.insert(mol->getAtomWithIdx(6));
    BOOST_TEST(!can_atoms_form_sgroup(atoms, *mol));
    // the atoms are now contiguous again
    atoms.insert(mol->getAtomWithIdx(5));
    BOOST_TEST(can_atoms_form_sgroup(atoms, *mol));
    // the atoms are still contiguous, but now there's only one connecting bond
    // since we've hit the end of the chain
    atoms.insert(mol->getAtomWithIdx(7));
    BOOST_TEST(!can_atoms_form_sgroup(atoms, *mol));

    // the empty set should always return false
    BOOST_TEST(!can_atoms_form_sgroup({}, *mol));

    // a single atom in the middle of the chain can form an S-group
    BOOST_TEST(can_atoms_form_sgroup({mol->getAtomWithIdx(3)}, *mol));
    // a single atom at the end of the chain can't since it only has one bond
    BOOST_TEST(!can_atoms_form_sgroup({mol->getAtomWithIdx(7)}, *mol));

    mol = rdkit_extensions::to_rdkit("CCC(CC)CC");
    // atom 1 has two bonds, so it can form an S-group
    BOOST_TEST(can_atoms_form_sgroup({mol->getAtomWithIdx(1)}, *mol));
    // but atom 2 has three bonds, so it can't
    BOOST_TEST(!can_atoms_form_sgroup({mol->getAtomWithIdx(2)}, *mol));

    // this set has three connecting bonds, so it can't form an S-group
    std::unordered_set<const RDKit::Atom*> atoms2(
        {mol->getAtomWithIdx(1), mol->getAtomWithIdx(2)});
    BOOST_TEST(!can_atoms_form_sgroup(atoms2, *mol));
    atoms2.insert(mol->getAtomWithIdx(3));
    BOOST_TEST(!can_atoms_form_sgroup(atoms2, *mol));
}

/**
 * Ensure that get_bonds_for_sgroup_atoms returns the expected bonds
 */
BOOST_AUTO_TEST_CASE(test_get_bonds_for_sgroup_atoms)
{
    auto mol = rdkit_extensions::to_rdkit("CCCCC");
    std::unordered_set<const RDKit::Atom*> atoms(
        {mol->getAtomWithIdx(2), mol->getAtomWithIdx(3)});
    auto bond_idxs = get_bonds_for_sgroup_atoms(atoms, *mol);
    BOOST_TEST(bond_idxs.size() == 2);
    std::unordered_set<unsigned int> bond_idxs_set(bond_idxs.begin(),
                                                   bond_idxs.end());
    BOOST_TEST(bond_idxs_set == std::unordered_set<unsigned int>({1, 3}));
}
/**
 * Ensure that get_bonds_within_sgroup returns the expected bonds
 */
BOOST_AUTO_TEST_CASE(test_get_bonds_within_sgroup)
{
    auto mol = rdkit_extensions::to_rdkit(S_GROUP_SMILES);
    auto s_group = RDKit::getSubstanceGroups(*mol)[0];
    auto bonds = get_bonds_within_sgroup(s_group);
    BOOST_TEST(bonds.size() == 3);
    std::unordered_set<unsigned int> bond_idxs;
    std::transform(bonds.begin(), bonds.end(),
                   std::inserter(bond_idxs, bond_idxs.begin()),
                   [](auto* bond) { return bond->getIdx(); });
    BOOST_TEST(bond_idxs == std::unordered_set<unsigned int>({3, 4, 5}));
}

} // namespace sketcher
} // namespace schrodinger
