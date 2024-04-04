/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: mol ops
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rdkit_extensions_molops

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <memory>
#include <vector>

#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/SubstanceGroup.h>

#include "schrodinger/rdkit_extensions/molops.h"

using schrodinger::rdkit_extensions::ExtractMolFragment;
using schrodinger::rdkit_extensions::removeHs;

BOOST_AUTO_TEST_CASE(test_partial_removeHs)
{
    auto ps = RDKit::SmilesParserParams();
    ps.removeHs = false;

    std::unique_ptr<RDKit::RWMol> mol{RDKit::SmilesToMol(
        "[H]C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]", ps)};

    // Remove 2 Hs on the starting CH3-. and all Hs on the other CH3-
    std::vector<unsigned> atoms_ids{0, 2, 10};

    removeHs(*mol, atoms_ids);

    BOOST_CHECK_EQUAL(RDKit::MolToSmiles(*mol), "[H]CC([H])([H])C([H])([H])C");

    for (auto atom : mol->atoms()) {
        BOOST_CHECK_EQUAL(atom->getIsotope(), 0);
    }
}

[[nodiscard]] static std::unique_ptr<RDKit::RWMol> get_test_mol()
{
    std::unique_ptr<RDKit::RWMol> mol{RDKit::SmilesToMol("CCCCCCCCCCCCCCC")};
    for (auto& atom : mol->atoms()) {
        atom->setProp("orig_idx", atom->getIdx());
    }

    for (auto& bond : mol->bonds()) {
        bond->setProp("orig_idx", bond->getIdx());
    }

    return mol;
}

namespace bdata = boost::unit_test::data;

BOOST_TEST_DONT_PRINT_LOG_VALUE(std::vector<unsigned int>)

BOOST_DATA_TEST_CASE(test_extract_atoms,
                     bdata::make(std::vector<std::vector<unsigned int>>{
                         {0, 2, 4, 6, 8, 10, 12},
                         {0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12},
                         {0, 2, 4, 6, 8, 10, 12, 100, 200, 300},
                     }),
                     selected_atoms)
{
    auto mol = get_test_mol();

    std::vector<unsigned int> expected_atoms{0, 2, 4, 6, 8, 10, 12};
    auto extracted_mol = ExtractMolFragment(*mol, selected_atoms);
    BOOST_TEST(extracted_mol->getNumAtoms() == expected_atoms.size());

    std::vector<unsigned int> extracted_atoms;
    for (auto& atom : extracted_mol->atoms()) {
        extracted_atoms.push_back(atom->template getProp<int>("orig_idx"));
    }

    BOOST_TEST(extracted_atoms == expected_atoms);
}

BOOST_AUTO_TEST_CASE(test_extract_bonds)
{
    auto mol = get_test_mol();

    for (auto& bond : mol->bonds()) {
        bond->setProp("test_prop", true);
    }

    for (auto& bond : mol->bonds()) {
        auto begin_idx = bond->getBeginAtomIdx();
        auto end_idx = bond->getEndAtomIdx();
        auto extracted_mol = ExtractMolFragment(*mol, {begin_idx, end_idx});

        BOOST_TEST(extracted_mol->getNumBonds() == 1);
        BOOST_TEST(extracted_mol->getBondWithIdx(0)->getProp<bool>(
                       "test_prop") == true);

        BOOST_TEST(extracted_mol->getNumAtoms() == 2);
        BOOST_TEST(extracted_mol->getAtomWithIdx(0)->getProp<int>("orig_idx") ==
                   begin_idx);
        BOOST_TEST(extracted_mol->getAtomWithIdx(1)->getProp<int>("orig_idx") ==
                   end_idx);
    }
}

struct selected_components {
    std::vector<bool> selected_atoms;
    std::vector<bool> selected_bonds;
};

[[nodiscard]] static selected_components
get_selected_components(::RDKit::RWMol& mol,
                        const std::vector<unsigned int>& atom_ids)
{
    const auto num_atoms = mol.getNumAtoms();
    std::vector<bool> selected_atoms(num_atoms);
    for (auto& atom_idx : atom_ids) {
        if (atom_idx < num_atoms) {
            selected_atoms[atom_idx] = true;
        }
    }

    std::vector<bool> selected_bonds(mol.getNumBonds());
    for (auto& bond : mol.bonds()) {
        if (selected_atoms[bond->getBeginAtomIdx()] &&
            selected_atoms[bond->getEndAtomIdx()]) {
            selected_bonds[bond->getIdx()] = true;
        }
    }

    return {std::move(selected_atoms), std::move(selected_bonds)};
}

BOOST_DATA_TEST_CASE(test_extract_substance_groups,
                     (bdata::make(std::vector<std::vector<unsigned int>>{
                          {0, 1, 2, 3, 4, 5},
                          {0, 2, 4, 6, 8, 10, 12},
                          {0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12},
                          {0, 2, 4, 6, 8, 10, 12, 100, 200, 300},
                      }) *
                      (bdata::make(std::vector<std::vector<unsigned int>>{
                           {},
                           {0, 1, 2, 3, 4},
                           {9, 10, 11},
                       }) *
                       bdata::make(std::vector<std::vector<unsigned int>>{
                           {},
                           {0, 1, 2},
                           {3, 4, 5},
                       }))) *
                         bdata::make(std::vector<std::vector<unsigned int>>{
                             {},
                             {3, 4},
                             {5, 6},
                         }),
                     selected_atoms_vec, sgroup_atoms, sgroup_bonds,
                     sgroup_patoms)
{
    auto mol = get_test_mol();

    ::RDKit::SubstanceGroup sgroup{mol.get(), "COP"};
    sgroup.setAtoms(sgroup_atoms);
    sgroup.setBonds(sgroup_bonds);
    sgroup.setParentAtoms(sgroup_patoms);
    ::RDKit::addSubstanceGroup(*mol, std::move(sgroup));

    auto has_selected_components = [&](auto& components, auto& ref_bitset) {
        return components.empty() ||
               std::all_of(
                   components.begin(), components.end(), [&](auto& idx) {
                       return idx < ref_bitset.size() && ref_bitset[idx];
                   });
    };

    auto extracted_mol = ExtractMolFragment(*mol, selected_atoms_vec);

    auto [selected_atoms, selected_bonds] =
        get_selected_components(*mol, selected_atoms_vec);
    // sgroup should only be copied if all components are selected
    auto flag = ::RDKit::getSubstanceGroups(*extracted_mol).size() == 1;
    BOOST_TEST(flag ==
               (has_selected_components(sgroup_atoms, selected_atoms) &&
                has_selected_components(sgroup_patoms, selected_atoms) &&
                has_selected_components(sgroup_bonds, selected_bonds)));

    // now make sure we copied the correct components
    if (flag) {
        auto& extracted_sgroup = ::RDKit::getSubstanceGroups(*extracted_mol)[0];
        for (auto& idx : extracted_sgroup.getAtoms()) {
            auto atom = extracted_mol->getAtomWithIdx(idx);
            BOOST_TEST(
                selected_atoms[atom->template getProp<int>("orig_idx")] ==
                true);
        }

        for (auto& idx : extracted_sgroup.getParentAtoms()) {
            auto atom = extracted_mol->getAtomWithIdx(idx);
            BOOST_TEST(
                selected_atoms[atom->template getProp<int>("orig_idx")] ==
                true);
        }

        for (auto& idx : extracted_sgroup.getBonds()) {
            auto bond = extracted_mol->getBondWithIdx(idx);
            BOOST_TEST(
                selected_bonds[bond->template getProp<int>("orig_idx")] ==
                true);
        }
    }
}

BOOST_DATA_TEST_CASE(test_extract_stereo_groups,
                     bdata::make(std::vector<std::vector<unsigned int>>{
                         {0, 1, 2, 3, 4, 5},
                         {0, 2, 4, 6, 8, 10, 12},
                         {0, 0, 2, 2, 4, 4, 6, 6, 8, 8, 10, 10, 12, 12},
                         {0, 2, 4, 6, 8, 10, 12, 100, 200, 300},
                     }) * (bdata::make(std::vector<std::vector<unsigned int>>{
                               {},
                               {0, 1, 2, 3, 4},
                               {9, 10, 11},
                           }) *
                           bdata::make(std::vector<std::vector<unsigned int>>{
                               {},
                               {0, 1, 2},
                               {3, 4, 5},
                           })),
                     selected_atoms_vec, stereo_group_atoms, stereo_group_bonds)
{
    auto mol = get_test_mol();

    std::vector<::RDKit::Atom*> sg_atoms;
    for (auto& idx : stereo_group_atoms) {
        sg_atoms.push_back(mol->getAtomWithIdx(idx));
    }
    std::vector<::RDKit::Bond*> sg_bonds;
    for (auto& idx : stereo_group_bonds) {
        sg_bonds.push_back(mol->getBondWithIdx(idx));
    }

    ::RDKit::StereoGroup stereo_group{::RDKit::StereoGroupType::STEREO_ABSOLUTE,
                                      std::move(sg_atoms), std::move(sg_bonds)};
    mol->setStereoGroups({std::move(stereo_group)});

    auto has_selected_components = [&](auto& components, auto& ref_bitset) {
        return components.empty() ||
               std::all_of(
                   components.begin(), components.end(), [&](auto& idx) {
                       return idx < ref_bitset.size() && ref_bitset[idx];
                   });
    };

    auto extracted_mol = ExtractMolFragment(*mol, selected_atoms_vec);

    auto [selected_atoms, selected_bonds] =
        get_selected_components(*mol, selected_atoms_vec);

    // stereo group should only be copied if all components are selected
    auto flag = extracted_mol->getStereoGroups().size() == 1;
    BOOST_TEST(flag ==
               (has_selected_components(stereo_group_atoms, selected_atoms) &&
                has_selected_components(stereo_group_bonds, selected_bonds)));

    // now make sure we copied the correct components
    if (flag) {
        auto& extracted_stereo_group = extracted_mol->getStereoGroups()[0];
        for (auto& atom : extracted_stereo_group.getAtoms()) {
            BOOST_TEST(
                selected_atoms[atom->template getProp<int>("orig_idx")] ==
                true);
        }

        for (auto& bond : extracted_stereo_group.getBonds()) {
            BOOST_TEST(
                selected_bonds[bond->template getProp<int>("orig_idx")] ==
                true);
        }
    }
}
