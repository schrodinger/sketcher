/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: coarse_grain
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rdkit_extensions_coarse_grain

#include <boost/test/unit_test.hpp>

#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <RDGeneral/RDLog.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>

#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"
#include "schrodinger/rdkit_extensions/coarse_grain.h"
#include "schrodinger/rdkit_extensions/cg_conversions.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/test/boost_checks.h"

using namespace schrodinger::rdkit_extensions;

BOOST_AUTO_TEST_CASE(TestBasicCoarseGrainMol)
{
    ::RDKit::RWMol cg_mol;
    auto monomer_idx1 = add_monomer(cg_mol, "A", 1, "PEPTIDE1");
    auto monomer_idx2 = add_monomer(cg_mol, "G");
    auto monomer_idx3 = add_monomer(cg_mol, "C");
    add_connection(cg_mol, monomer_idx1, monomer_idx2);
    add_connection(cg_mol, monomer_idx2, monomer_idx3);
    assign_chains(cg_mol);

    BOOST_CHECK_EQUAL(to_string(cg_mol, Format::HELM),
                      "PEPTIDE1{A.G.C}$$$$V2.0");
    BOOST_CHECK_EQUAL(to_string(cg_mol, Format::FASTA), ">\nAGC\n");

    // making it a cyclic peptide
    add_connection(cg_mol, 2, 0);
    assign_chains(cg_mol);
    BOOST_CHECK_EQUAL(to_string(cg_mol, Format::HELM),
                      "PEPTIDE1{A.G.C}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$V2.0");
    // This is a cyclic peptide, should FASTA throw?
    // BOOST_CHECK_THROW(to_string(cg_mol, Format::FASTA),
    // std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestBranchesCoarseGrain)
{
    ::RDKit::RWMol cg_mol;
    auto monomer_idx1 = add_monomer(cg_mol, "A", 1, "PEPTIDE1");
    auto monomer_idx2 = add_monomer(cg_mol, "G");
    auto monomer_branch = add_monomer(cg_mol, "A");
    auto monomer_idx3 = add_monomer(cg_mol, "C");
    auto monomer_idx4 = add_monomer(cg_mol, "T");

    add_connection(cg_mol, monomer_idx1, monomer_idx2);
    add_connection(cg_mol, monomer_idx2, monomer_idx3);
    add_connection(cg_mol, monomer_idx3, monomer_idx4);
    add_connection(cg_mol, monomer_idx2, monomer_branch,
                   ConnectionType::SIDECHAIN);
    assign_chains(cg_mol);
    BOOST_CHECK_EQUAL(to_string(cg_mol, Format::HELM),
                      "PEPTIDE1{A.G(A)C.T}$$$$V2.0");

    // a monomer not in connectivity order will be placed
    // incorrectly in the HELM string, so for now we throw
    auto monomer_branch2 = add_monomer(cg_mol, "A");
    add_connection(cg_mol, monomer_idx3, monomer_branch2,
                   ConnectionType::SIDECHAIN);
    assign_chains(cg_mol);
    BOOST_CHECK_EQUAL(to_string(cg_mol, Format::HELM),
                      "PEPTIDE1{A.G(A)C.T(A)}$$$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestMultipleChainsCoarseGrainMol)
{
    ::RDKit::RWMol cg_mol;
    auto monomer_idx1 = add_monomer(cg_mol, "A", 1, "PEPTIDE1");
    auto monomer_idx2 = add_monomer(cg_mol, "G");
    auto monomer_idx3 = add_monomer(cg_mol, "C");
    add_connection(cg_mol, monomer_idx1, monomer_idx2);
    add_connection(cg_mol, monomer_idx2, monomer_idx3);

    auto monomer_idx4 = add_monomer(cg_mol, "T", 1, "PEPTIDE2");
    auto monomer_idx5 = add_monomer(cg_mol, "C");
    auto monomer_idx6 = add_monomer(cg_mol, "A");
    add_connection(cg_mol, monomer_idx4, monomer_idx5);
    add_connection(cg_mol, monomer_idx5, monomer_idx6);

    add_connection(cg_mol, monomer_idx3, monomer_idx4,
                   ConnectionType::SIDECHAIN);

    assign_chains(cg_mol);
    BOOST_CHECK_EQUAL(
        to_string(cg_mol, Format::HELM),
        "PEPTIDE1{A.G.C}|PEPTIDE2{T.C.A}$PEPTIDE1,PEPTIDE2,3:R3-1:R1$$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestAtomisticSmilesToCGString)
{
    {
        std::unique_ptr<::RDKit::RWMol> mol(::RDKit::SmilesToMol(
            "CC[C@H](C)[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@"
            "H](CCC(=O)O)NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CO)NC(=O)["
            "C@@H](NC(=O)[C@H](Cc1ccc(O)cc1)NC(C)=O)[C@@H](C)O)[C@@H](C)CC)C(="
            "O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@H](C(=O)N[C@@H](CCC(="
            "O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC("
            "N)=O)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(N)=O)C("
            "=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C("
            "=O)N[C@@H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@"
            "H](CC(C)C)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H]"
            "(Cc1c[nH]c2ccccc12)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H]("
            "CCCCN)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(N)=O)C(=O)N["
            "C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccccc1)C(N)=O)[C@@H](C)"
            "O"));
        auto cg_mol = atomistic_to_cg(*mol);
        BOOST_CHECK_EQUAL(to_string(*cg_mol, Format::HELM),
                          "PEPTIDE1{I.L.E.E.I.L.S.T.Y.K.K.T.E.E.Q.Q.K.K.N.E.E."
                          "E.L.K.K.L.E.E.W.A.K.K.W.N.W.F}$$$$V2.0");
        BOOST_CHECK_EQUAL(to_string(*cg_mol, Format::FASTA),
                          ">\nILEEILSTYKKTEEQQKKNEEELKKLEEWAKKWNWF\n");
    }

    {
        // cyclic peptide with non-standard (or just unrecognized) monomer
        std::unique_ptr<::RDKit::RWMol> mol(::RDKit::SmilesToMol(
            "CC(C)C[C@@H]1NC(=O)[C@H](CCCN)NC(=O)[C@H](C(C)C)NC(=O)[C@@H]"
            "2CCCN2C(=O)[C@@H](Cc2ccccc2)NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCCN)"
            "NC(=O)[C@H](C(C)C)NC(=O)[C@@H]2CCCN2C(=O)[C@@H](Cc2ccccc2)NC1=O"));
        auto cg_mol = atomistic_to_cg(*mol);
        BOOST_CHECK_EQUAL(to_string(*cg_mol, Format::HELM),
                          "PEPTIDE1{L.[NCCCC(N)C(N)=O].V.P.F.L.[NCCCC(N)C(N)=O]"
                          ".V.P.F}$PEPTIDE1,PEPTIDE1,1:R2-10:R1$$$V2.0");
        // SMILES monomer not allowed in FASTA
        BOOST_CHECK_THROW(to_string(*cg_mol, Format::FASTA),
                          std::invalid_argument);
    }
}
