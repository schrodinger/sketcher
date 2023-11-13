/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: mol ops
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rdkit_extensions_convert

#include <boost/test/unit_test.hpp>

#include <memory>
#include <vector>

#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

#include "schrodinger/rdkit_extensions/molops.h"

using namespace schrodinger;

BOOST_AUTO_TEST_CASE(test_partial_removeHs)
{
    auto ps = RDKit::SmilesParserParams();
    ps.removeHs = false;

    std::unique_ptr<RDKit::RWMol> mol{RDKit::SmilesToMol(
        "[H]C([H])([H])C([H])([H])C([H])([H])C([H])([H])[H]", ps)};

    // Remove 2 Hs on the starting CH3-. and all Hs on the other CH3-
    std::vector<unsigned> atoms_ids{0, 2, 10};

    rdkit_extensions::removeHs(*mol, atoms_ids);

    BOOST_CHECK_EQUAL(RDKit::MolToSmiles(*mol), "[H]CC([H])([H])C([H])([H])C");

    for (auto atom : mol->atoms()) {
        BOOST_CHECK_EQUAL(atom->getIsotope(), 0);
    }
}
