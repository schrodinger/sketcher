#define BOOST_TEST_MODULE Test_Sketcher

#include <algorithm>
#include <vector>

#include <boost/test/data/test_case.hpp>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/monomer_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

namespace bdata = boost::unit_test::data;

const std::vector<std::string> HELM_IN = {
    // two neighboring cysteines, but no disulfidea
    "PEPTIDE1{C.C}$$$$V2.0",
    // two neighboring cysteines with a disulfide
    "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0",
    // a disulfide, but between two non-neighboring cysteines
    "PEPTIDE1{C.A.C}$PEPTIDE1,PEPTIDE1,1:R3-3:R3$$$V2.0"};
const std::vector<int> EXP_NUM_ATOMS = {2, 2, 3};
const std::vector<int> EXP_NUM_CONNECTIONS = {1, 1, 3};
const std::vector<int> EXP_NUM_SECONDARY_CONNECTIONS = {0, 1, 0};

/**
 * Make sure that create_graphics_items_for_mol() creates the expected number of
 * graphics items for monomeric models, including returning two separate
 * graphics items for bonds with a secondary connections. Also ensure that
 * get_model_objects_for_graphics_items() can return the corresponding model
 * objects.
 */
BOOST_DATA_TEST_CASE(test_create_graphics_items_for_mol_monomeric,
                     bdata::make(HELM_IN) ^ bdata::make(EXP_NUM_ATOMS) ^
                         bdata::make(EXP_NUM_CONNECTIONS) ^
                         bdata::make(EXP_NUM_SECONDARY_CONNECTIONS),
                     helm_in, exp_num_atoms, exp_num_connections,
                     exp_num_secondary_connections)
{
    AtomDisplaySettings atom_display_settings;
    BondDisplaySettings bond_display_settings;
    Fonts fonts;

    auto mol = rdkit_extensions::to_rdkit(helm_in);
    prepare_mol(*mol);
    for (auto* atom : mol->atoms()) {
        set_atom_monomeric(atom);
    }
    auto [all_items, atom_items, bond_items, secondary_connection_items,
          sgroup_items] =
        create_graphics_items_for_mol(mol.get(), fonts, atom_display_settings,
                                      bond_display_settings);
    BOOST_TEST(all_items.size() == exp_num_atoms + exp_num_connections +
                                       exp_num_secondary_connections);
    BOOST_TEST(atom_items.size() == exp_num_atoms);
    BOOST_TEST(bond_items.size() == exp_num_connections);
    BOOST_TEST(secondary_connection_items.size() ==
               exp_num_secondary_connections);
    BOOST_TEST(sgroup_items.empty());

    std::unordered_set<QGraphicsItem*> all_items_set;
    std::copy(all_items.begin(), all_items.end(),
              std::inserter(all_items_set, all_items_set.end()));
    auto [atoms, bonds, secondary_connections, sgroups, non_molecular_objects] =
        get_model_objects_for_graphics_items(all_items_set);
    BOOST_TEST(atoms.size() == exp_num_atoms);
    BOOST_TEST(bonds.size() == exp_num_connections);
    BOOST_TEST(secondary_connections.size() == exp_num_secondary_connections);
    BOOST_TEST(sgroups.empty());
    BOOST_TEST(non_molecular_objects.empty());

    // we don't have a Scene in the test, so we have to delete the graphics
    // items manually
    for (auto item : all_items) {
        delete item;
    }
}

std::vector<int> get_linear_custom_connector_lanes(const std::string& helm)
{
    auto mol = rdkit_extensions::to_rdkit(helm);
    rdkit_extensions::compute_monomer_mol_coords(
        *mol, /* render_peptides_linearly = */ true);
    prepare_mol(*mol);
    for (auto* atom : mol->atoms()) {
        set_atom_monomeric(atom);
    }

    AtomDisplaySettings atom_display_settings;
    BondDisplaySettings bond_display_settings;
    Fonts fonts;
    auto [all_items, atom_items, bond_items, secondary_connection_items,
          sgroup_items] =
        create_graphics_items_for_mol(mol.get(), fonts, atom_display_settings,
                                      bond_display_settings,
                                      /* is_dark_mode = */ false,
                                      /* draw_attachment_points = */ true,
                                      /* render_peptides_linearly = */ true);

    std::vector<const RDKit::Bond*> custom_bonds;
    for (const auto* bond : mol->bonds()) {
        if (bond->hasProp(CUSTOM_BOND)) {
            custom_bonds.push_back(bond);
        }
    }
    std::sort(custom_bonds.begin(), custom_bonds.end(),
              [](const auto* lhs, const auto* rhs) {
                  const auto lhs_start =
                      std::min(lhs->getBeginAtomIdx(), lhs->getEndAtomIdx());
                  const auto rhs_start =
                      std::min(rhs->getBeginAtomIdx(), rhs->getEndAtomIdx());
                  return lhs_start < rhs_start;
              });
    std::vector<int> lanes;
    for (const auto* bond : custom_bonds) {
        const auto* connector =
            dynamic_cast<const MonomerConnectorItem*>(bond_items.at(bond));
        BOOST_REQUIRE(connector != nullptr);
        lanes.push_back(connector->getConnectorLane());
    }

    for (auto* item : all_items) {
        delete item;
    }
    return lanes;
}

BOOST_AUTO_TEST_CASE(test_interlocking_linear_peptide_connectors_use_both_sides)
{
    const auto lanes = get_linear_custom_connector_lanes(
        "PEPTIDE1{C.K.F.K.L.C.C.T.Y.A.G.C}$PEPTIDE1,PEPTIDE1,"
        "1:R3-7:R3|PEPTIDE1,PEPTIDE1,6:R3-12:R3$$$V2.0");
    BOOST_REQUIRE_EQUAL(lanes.size(), 2);
    BOOST_CHECK_EQUAL(lanes[0], 1);
    BOOST_CHECK_EQUAL(lanes[1], -1);
}

BOOST_AUTO_TEST_CASE(test_nested_linear_peptide_connectors_use_both_sides)
{
    const auto lanes = get_linear_custom_connector_lanes(
        "PEPTIDE1{[Glp].G.V.C.C.G.Y.K.L.C.H.[dD].C}$PEPTIDE1,PEPTIDE1,"
        "5:R3-10:R3|PEPTIDE1,PEPTIDE1,4:R3-13:R3$$$");
    BOOST_REQUIRE_EQUAL(lanes.size(), 2);
    // The shorter inner bond is assigned first (top), so the outer bond is
    // routed to the next available lane (bottom).
    BOOST_CHECK_EQUAL(lanes[0], -1);
    BOOST_CHECK_EQUAL(lanes[1], 1);
}

BOOST_AUTO_TEST_CASE(test_overlapping_linear_peptide_connectors_use_more_lanes)
{
    const auto lanes = get_linear_custom_connector_lanes(
        "PEPTIDE1{C.C.[dN].C.S.A.K.W.C.R.D.H.S.R.C.C.[am]}$"
        "PEPTIDE1,PEPTIDE1,1:R3-9:R3|PEPTIDE1,PEPTIDE1,2:R3-15:R3|"
        "PEPTIDE1,PEPTIDE1,4:R3-16:R3$$$");
    BOOST_REQUIRE_EQUAL(lanes.size(), 3);
    BOOST_CHECK_EQUAL(lanes[0], 1);
    BOOST_CHECK_EQUAL(lanes[1], 2);
    BOOST_CHECK_EQUAL(lanes[2], -1);
}

BOOST_AUTO_TEST_CASE(test_linear_peptide_r2_r1_connection_uses_a_lane)
{
    const auto lanes = get_linear_custom_connector_lanes(
        "PEPTIDE1{A.A.G.F.P.V.F.F}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$");
    BOOST_REQUIRE_EQUAL(lanes.size(), 1);
    BOOST_CHECK_EQUAL(lanes[0], 1);
}

} // namespace sketcher
} // namespace schrodinger
