/* -------------------------------------------------------------------------
 * Tests demonstrating RDKit behaviors with wiggly bonds (unknown
 * stereochemistry). These tests document the specific RDKit behaviors that
 * necessitate wiggly bond preservation code in the sketcher.
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_MODULE rdkit_extensions_wiggly_bond_behaviors

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

/**
 * Demonstrate that CXSMILES wiggly bonds lack _MolFileBondCfg property
 * (which prepare_mol checks to identify wiggly bonds for preservation).
 * This test also shows that insertMol preserves BondDir but not the property.
 */
BOOST_AUTO_TEST_CASE(test_cxsmiles_wiggly_bonds_lack_molfile_property)
{
    // Create a molecule with a wiggly bond using CXSMILES
    const std::string cxsmiles = "CCC(C)N |w:2.3|";
    RDKit::SmilesParserParams params;
    params.sanitize = true;
    params.removeHs = false;
    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol(cxsmiles, params));
    BOOST_REQUIRE(mol != nullptr);

    // Verify CXSMILES roundtrip preserves wiggly bond
    std::string cxsmiles_output = RDKit::MolToCXSmiles(*mol);
    BOOST_TEST(cxsmiles_output == cxsmiles);

    // Get the wiggly bond (bond index 2: between atom 2 and atom 3)
    auto* bond = mol->getBondWithIdx(2);
    BOOST_REQUIRE(bond != nullptr);

    // The key issue: CXSMILES input does NOT set _MolFileBondCfg property
    // (unlike MDL input which sets both BondDir AND the property)
    int cfg{-1};
    bond->getPropIfPresent(RDKit::common_properties::_MolFileBondCfg, cfg);
    BOOST_TEST(cfg == -1); // property not set for CXSMILES input

    // RDKit's CXSMILES parser also doesn't set BondDir to UNKNOWN initially.
    // The wiggly bond information is only preserved in the CXSMILES extension.
    // When coordinates are generated and the molecule is processed, BondDir
    // gets set.

    // However, we can manually set BondDir to UNKNOWN (like mol_model.cpp does)
    // to simulate what happens when the mol is copied
    bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
    BOOST_TEST(bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);

    // Test insertMol preserves BondDir
    RDKit::RWMol target_mol;
    target_mol.insertMol(*mol);

    auto* dest_bond = target_mol.getBondWithIdx(2);
    BOOST_REQUIRE(dest_bond != nullptr);
    BOOST_TEST(dest_bond->getBondDir() == RDKit::Bond::BondDir::UNKNOWN);

    // But the property still isn't there - this is why we need to set it in
    // mol_model.cpp
    cfg = -1;
    dest_bond->getPropIfPresent(RDKit::common_properties::_MolFileBondCfg, cfg);
    BOOST_TEST(cfg == -1);
}
