/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: coarse_grain
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_MODULE rdkit_extensions_coarse_grain

#include <boost/test/unit_test.hpp>

#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/RDGeneral/RDLog.h>

#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/rdkit_extensions/atomistic_conversions.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "test_common.h"

using namespace schrodinger::rdkit_extensions;

BOOST_AUTO_TEST_CASE(TestBasicCoarseGrainMol)
{
    ::RDKit::RWMol monomer_mol;
    auto monomer_idx1 = addMonomer(monomer_mol, "A", 1, "PEPTIDE1");
    auto monomer_idx2 = addMonomer(monomer_mol, "G");
    auto monomer_idx3 = addMonomer(monomer_mol, "C");
    addConnection(monomer_mol, monomer_idx1, monomer_idx2);
    addConnection(monomer_mol, monomer_idx2, monomer_idx3);
    assignChains(monomer_mol);

    BOOST_CHECK_EQUAL(to_string(monomer_mol, Format::HELM),
                      "PEPTIDE1{A.G.C}$$$$V2.0");
    BOOST_CHECK_EQUAL(to_string(monomer_mol, Format::FASTA), ">\nAGC\n");

    // making it a cyclic peptide
    addConnection(monomer_mol, 2, 0);
    assignChains(monomer_mol);
    BOOST_CHECK_EQUAL(to_string(monomer_mol, Format::HELM),
                      "PEPTIDE1{A.G.C}$PEPTIDE1,PEPTIDE1,3:R2-1:R1$$$V2.0");
    // This is a cyclic peptide, should FASTA throw?
    // BOOST_CHECK_THROW(to_string(monomer_mol, Format::FASTA),
    // std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestBranchesCoarseGrain)
{
    ::RDKit::RWMol monomer_mol;
    auto monomer_idx1 = addMonomer(monomer_mol, "A", 1, "PEPTIDE1");
    auto monomer_idx2 = addMonomer(monomer_mol, "C");
    auto monomer_branch = addMonomer(monomer_mol, "A");
    auto monomer_idx3 = addMonomer(monomer_mol, "C");
    auto monomer_idx4 = addMonomer(monomer_mol, "T");

    addConnection(monomer_mol, monomer_idx1, monomer_idx2);
    addConnection(monomer_mol, monomer_idx2, monomer_idx3);
    addConnection(monomer_mol, monomer_idx3, monomer_idx4);
    addConnection(monomer_mol, monomer_idx2, monomer_branch,
                  ConnectionType::SIDECHAIN);
    assignChains(monomer_mol);
    BOOST_CHECK_EQUAL(to_string(monomer_mol, Format::HELM),
                      "PEPTIDE1{A.C(A)C.T}$$$$V2.0");

    auto monomer_branch2 = addMonomer(monomer_mol, "A");
    addConnection(monomer_mol, monomer_idx3, monomer_branch2,
                  ConnectionType::SIDECHAIN);
    assignChains(monomer_mol);
    BOOST_CHECK_EQUAL(to_string(monomer_mol, Format::HELM),
                      "PEPTIDE1{A.C(A)C(A)T}$$$$V2.0");
    // TODO (SHARED-10771): Fix HELM writer, there should be a `.` bewteen
    // consecutive branches
}

BOOST_AUTO_TEST_CASE(TestMultipleChainsCoarseGrainMol)
{
    ::RDKit::RWMol monomer_mol;
    auto monomer_idx1 = addMonomer(monomer_mol, "A", 1, "PEPTIDE1");
    auto monomer_idx2 = addMonomer(monomer_mol, "G");
    auto monomer_idx3 = addMonomer(monomer_mol, "C");
    addConnection(monomer_mol, monomer_idx1, monomer_idx2);
    addConnection(monomer_mol, monomer_idx2, monomer_idx3);

    auto monomer_idx4 = addMonomer(monomer_mol, "T", 1, "PEPTIDE2");
    auto monomer_idx5 = addMonomer(monomer_mol, "C");
    auto monomer_idx6 = addMonomer(monomer_mol, "A");
    addConnection(monomer_mol, monomer_idx4, monomer_idx5);
    addConnection(monomer_mol, monomer_idx5, monomer_idx6);

    addConnection(monomer_mol, monomer_idx3, monomer_idx4,
                  ConnectionType::SIDECHAIN);
    assignChains(monomer_mol);

    BOOST_CHECK(monomer_mol.getAtomWithIdx(monomer_idx4)
                    ->getProp<bool>(BRANCH_MONOMER));
    BOOST_CHECK(!monomer_mol.getAtomWithIdx(monomer_idx3)
                     ->getProp<bool>(BRANCH_MONOMER));
    BOOST_CHECK_EQUAL(
        to_string(monomer_mol, Format::HELM),
        "PEPTIDE1{A.G.C}|PEPTIDE2{T.C.A}$PEPTIDE1,PEPTIDE2,3:R3-1:R1$$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestAtomisticSmilesToCGString)
{
    {
        std::unique_ptr<::RDKit::RWMol> mol(::RDKit::SmilesToMol(
            "CC[C@H](C)[C@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@"
            "H](CCC(=O)O)NC(=O)[C@@H](NC(=O)[C@H](CC(C)C)NC(=O)[C@H](CO)NC(=O)["
            "C@@H](NC(=O)[C@@H](N)Cc1ccc(O)cc1)[C@@H](C)O)[C@@H](C)CC)C(=O)N[C@"
            "@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)O)C(="
            "O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(N)=O)C(=O)N[C@@H](CCC(N)=O)C(="
            "O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@"
            "H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@"
            "H](CC(C)C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CC(C)"
            "C)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](CCC(=O)O)C(=O)N[C@@H](Cc1c["
            "nH]c2ccccc12)C(=O)N[C@@H](C)C(=O)N[C@@H](CCCCN)C(=O)N[C@@H](CCCCN)"
            "C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CC(N)=O)C(=O)N[C@@H]("
            "Cc1c[nH]c2ccccc12)C(=O)N[C@@H](Cc1ccccc1)C(=O)O)[C@@H](C)O"));
        bool use_residue_info = false;
        auto monomer_mol = toMonomeric(*mol, use_residue_info);
        BOOST_CHECK_EQUAL(to_string(*monomer_mol, Format::HELM),
                          "PEPTIDE1{Y.T.S.L.I.E.E.L.I.K.K.T.E.E.Q.Q.K.K.N.E.E."
                          "E.L.K.K.L.E.E.W.A.K.K.W.N.W.F}$$$$V2.0");
        BOOST_CHECK_EQUAL(to_string(*monomer_mol, Format::FASTA),
                          ">\nYTSLIEELIKKTEEQQKKNEEELKKLEEWAKKWNWF\n");
    }

    {
        // cyclic peptide
        std::unique_ptr<::RDKit::RWMol> mol(::RDKit::SmilesToMol(
            "CCC(C)[C@@H]1NC(=O)C2CCCN2C(=O)[C@H](CCCCN)NC(=O)[C@H](CC(N)=O)NC("
            "=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](CC(C)C)"
            "NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H]2CCCN2C1=O"));
        bool use_residue_info = false;
        auto monomer_mol = toMonomeric(*mol, use_residue_info);
        BOOST_CHECK_EQUAL(to_string(*monomer_mol, Format::HELM),
                          "PEPTIDE1{P.P.F.L.W.L.N.K.P.I}$PEPTIDE1,PEPTIDE1,10:"
                          "R2-1:R1$$$V2.0");
    }
}

BOOST_AUTO_TEST_CASE(Test_annotated_toMonomeric)
{
    auto atomistic_mol =
        RDKit::v2::FileParsers::MolFromPDBFile(testfile_path("1dng.pdb"));
    auto monomer_mol = toMonomeric(*atomistic_mol);

    BOOST_CHECK_EQUAL(to_string(*monomer_mol, Format::HELM),
                      "PEPTIDE1{Q.A.P.A.Y.E.E.A.A.E.E.L.A.K.S}$$$$V2.0");
}

BOOST_AUTO_TEST_CASE(Test_toMonomeric)
{
    // More tests validation tests are available in python, these are to ensure
    // memtest is ran
    {
        auto monomer_mol =
            to_rdkit("PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0",
                     Format::HELM);
        auto atomistic_mol = toAtomistic(*monomer_mol);
        BOOST_CHECK_EQUAL(::RDKit::MolToSmiles(*atomistic_mol),
                          "C[C@@H]1NC(=O)[C@H](C)NC(=O)[C@@H](N)CSSC[C@@H](C(="
                          "O)O)NC(=O)[C@H](C)NC1=O");
    }
    {
        auto monomer_mol = to_rdkit("PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,"
                                    "PEPTIDE2,2:R3-1:R1$$$V2.0");
        auto atomistic_mol = toAtomistic(*monomer_mol);
        BOOST_CHECK_EQUAL(
            ::RDKit::MolToSmiles(*atomistic_mol),
            "N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)N[C@@H](CS)C(=O)N[C@@H](CCC(="
            "O)O)C(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)NCC(=O)O");
    }
    {
        // error condition -- invalid monomer
        auto monomer_mol = to_rdkit("PEPTIDE1{F.Y.Z.G.R.L}$$$$V2.0");
        BOOST_CHECK_THROW(toAtomistic(*monomer_mol), std::out_of_range);
    }
    {
        // error condition -- A does not have an R3 attachment point
        auto monomer_mol =
            to_rdkit("PEPTIDE1{A.A.A.A.A}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0");
        BOOST_CHECK_THROW(toAtomistic(*monomer_mol), std::runtime_error);
    }
    {
        // error condition -- A/C/G/T are the nitrogenous base thus only have 1
        // attachment point and cannot be in the backbone
        auto monomer_mol = to_rdkit("RNA1{A.C.G.T}$$$$V2.0");
        BOOST_CHECK_THROW(toAtomistic(*monomer_mol), std::runtime_error);
    }
}

BOOST_AUTO_TEST_CASE(Test_reordering_residues)
{
    ::RDKit::RWMol monomer_mol;
    addMonomer(monomer_mol, "D", 3, "PEPTIDE1");
    addMonomer(monomer_mol, "C", 1, "PEPTIDE1");
    addMonomer(monomer_mol, "B", 4, "PEPTIDE1");
    addMonomer(monomer_mol, "A", 2, "PEPTIDE1");
    addConnection(monomer_mol, 3, 2);
    addConnection(monomer_mol, 2, 1);
    addConnection(monomer_mol, 1, 0);
    assignChains(monomer_mol);

    BOOST_CHECK_EQUAL(to_string(monomer_mol, Format::HELM),
                      "PEPTIDE1{A.B.C.D}$$$$V2.0");
}