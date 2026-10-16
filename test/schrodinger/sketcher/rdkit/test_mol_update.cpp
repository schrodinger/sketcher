
#define BOOST_TEST_MODULE mol_update

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/constants.h"
#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>

namespace schrodinger
{
namespace sketcher
{

namespace bdata = boost::unit_test::data;
namespace tt = boost::test_tools;

std::string LIST_QUERY = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 3 2 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -14.969697 7.242424 0.000000 0
M  V30 2 [C,N] -13.670659 7.992424 0.000000 0
M  V30 3 C -12.371621 7.242424 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$)";

// SKETCH-2710: a molecule with query bonds (V3000 bond type 8 = "any") and a
// stereocenter must not throw when update_molecule_on_change is called. The
// CIPLabeler cannot handle non-integer (UNSPECIFIED) bond orders produced by
// query bonds, so CIP label assignment must be skipped for such molecules.
std::string QUERY_BOND_WITH_STEREO = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 11 11 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C 1.275976 0.000000 0.000000 0
M  V30 2 C 0.394298 1.213525 0.000000 0
M  V30 3 C -1.032286 0.750000 0.000000 0
M  V30 4 C -1.032286 -0.750000 0.000000 0
M  V30 5 C 0.394298 -1.213525 0.000000 0
M  V30 6 C -2.245812 1.631678 0.000000 0
M  V30 7 C -2.089019 3.123461 0.000000 0
M  V30 8 C -0.718701 3.733566 0.000000 0
M  V30 9 C -3.302545 4.005139 0.000000 0
M  V30 10 N 0.494824 2.851888 0.000000 0
M  V30 11 C -0.561908 5.225349 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 8 1 2
M  V30 2 8 2 3
M  V30 3 8 3 4
M  V30 4 8 4 5
M  V30 5 8 5 1
M  V30 6 1 3 6
M  V30 7 2 6 7
M  V30 8 1 7 8
M  V30 9 1 7 9
M  V30 10 1 8 10
M  V30 11 1 8 11 CFG=1
M  V30 END BOND
M  V30 BEGIN COLLECTION
M  V30 MDLV30/STEABS ATOMS=(1 8)
M  V30 END COLLECTION
M  V30 END CTAB
M  END
)";

BOOST_AUTO_TEST_CASE(test_update_mol_query_bonds_with_stereo_no_throw)
{
    auto mol = rdkit_extensions::to_rdkit(
        QUERY_BOND_WITH_STEREO, rdkit_extensions::Format::MDL_MOLV3000);
    BOOST_REQUIRE(mol != nullptr);
    // prepare_mol applies wedge bond directions from V3000 CFG fields, which
    // sets CHI tags on atom 3 — needed to trigger CIPLabeler traversal.
    prepare_mol(*mol);
    BOOST_REQUIRE_NO_THROW(update_molecule_on_change(*mol));
}

/**
 * Make sure that prepare_mol flags list queries as dummy atoms
 */
BOOST_AUTO_TEST_CASE(test_prepare_mol_list_query)
{
    auto mol = rdkit_extensions::to_rdkit(
        LIST_QUERY, rdkit_extensions::Format::MDL_MOLV3000);
    BOOST_TEST(mol->getNumAtoms() == 3);

    prepare_mol(*mol);

    // the second atom is a list query, check that it's a dummy atom
    BOOST_TEST(mol->getAtomWithIdx(1)->getAtomicNum() ==
               rdkit_extensions::DUMMY_ATOMIC_NUMBER);
}

static std::vector<RDKit::Bond::BondDir> get_bond_dirs(const RDKit::ROMol& mol)
{
    std::vector<RDKit::Bond::BondDir> dirs;
    for (auto bond : mol.bonds()) {
        dirs.push_back(bond->getBondDir());
    }
    return dirs;
}

/**
 * SKETCH-2872: prepare_mol should keep wedges already present on an input mol
 * with coordinates, and only wedge chiral centers that lack one
 */
BOOST_AUTO_TEST_CASE(test_prepare_mol_keeps_input_wedges)
{
    const std::string smiles = "C[C@H](O)C[C@H](N)C";
    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(smiles));
    RDDepict::compute2DCoords(*mol);

    // Default wedging for these coordinates
    RDKit::RWMol default_mol(*mol);
    prepare_mol(default_mol);
    auto default_dirs = get_bond_dirs(default_mol);

    // Manually wedge a bond on atom 1 that isn't wedged by default (bonds 1
    // and 2 both begin at atom 1), and leave atom 4 unwedged
    unsigned int wedged_idx =
        default_dirs[1] == RDKit::Bond::BondDir::NONE ? 1 : 2;
    BOOST_REQUIRE(default_dirs[wedged_idx] == RDKit::Bond::BondDir::NONE);
    RDKit::Chirality::wedgeBond(mol->getBondWithIdx(wedged_idx), 1,
                                &mol->getConformer());
    auto wedged_dir = mol->getBondWithIdx(wedged_idx)->getBondDir();
    BOOST_REQUIRE(wedged_dir != RDKit::Bond::BondDir::NONE);

    prepare_mol(*mol);
    for (auto bond : mol->bonds()) {
        auto idx = bond->getIdx();
        if (idx == wedged_idx) {
            // the manual wedge is kept
            BOOST_TEST(bond->getBondDir() == wedged_dir);
        } else if (bond->getBeginAtomIdx() == 1 || bond->getEndAtomIdx() == 1) {
            // atom 1 doesn't get a second wedge
            BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::NONE);
        } else {
            // atom 4 is wedged as normal
            BOOST_TEST(bond->getBondDir() == default_dirs[idx]);
        }
    }

    // If coordinates have to be generated, the input wedges are meaningless and
    // should be recalculated from scratch
    RDKit::RWMol no_coords_mol(*mol);
    no_coords_mol.clearConformers();
    std::unique_ptr<RDKit::RWMol> expected_mol(RDKit::SmilesToMol(smiles));
    prepare_mol(no_coords_mol);
    prepare_mol(*expected_mol);
    BOOST_TEST(get_bond_dirs(no_coords_mol) == get_bond_dirs(*expected_mol),
               tt::per_element());
}

} // namespace sketcher
} // namespace schrodinger