/* -------------------------------------------------------------------------
 * Tests class schrodinger::rdkit_extensions:: coarse_grain
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE rdkit_extensions_coarse_grain

#include <boost/test/unit_test.hpp>

#include <GraphMol/FileParsers/FileParsers.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>
#include <GraphMol/SmilesParse/SmilesWrite.h>
#include <RDGeneral/RDLog.h>

#include "schrodinger/rdkit_extensions/capture_rdkit_log.h"
#include "schrodinger/rdkit_extensions/coarse_grain.h"
#include "schrodinger/rdkit_extensions/cg_conversions.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/test/boost_checks.h"
#include "schrodinger/test/testfiles.h"

using namespace schrodinger::rdkit_extensions;
using schrodinger::test::mmshare_testfile;

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
    auto monomer_idx2 = add_monomer(cg_mol, "C");
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
                      "PEPTIDE1{A.C(A)C.T}$$$$V2.0");

    auto monomer_branch2 = add_monomer(cg_mol, "A");
    add_connection(cg_mol, monomer_idx3, monomer_branch2,
                   ConnectionType::SIDECHAIN);
    assign_chains(cg_mol);
    BOOST_CHECK_EQUAL(to_string(cg_mol, Format::HELM),
                      "PEPTIDE1{A.C(A)C(A)T}$$$$V2.0");
    // TODO (SHARED-10771): Fix HELM writer, there should be a `.` bewteen
    // consecutive branches
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

    BOOST_CHECK(
        cg_mol.getAtomWithIdx(monomer_idx4)->getProp<bool>(BRANCH_MONOMER));
    BOOST_CHECK(
        !cg_mol.getAtomWithIdx(monomer_idx3)->getProp<bool>(BRANCH_MONOMER));
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
        bool use_residue_info = false;
        auto cg_mol = atomistic_to_cg(*mol, use_residue_info);
        BOOST_CHECK_EQUAL(to_string(*cg_mol, Format::HELM),
                          "PEPTIDE1{Y.T.S.L.I.E.E.L.I.K.K.T.E.E.Q.Q.K.K.N.E.E."
                          "E.L.K.K.L.E.E.W.A.K.K.W.N.W.F}$$$$V2.0");
        BOOST_CHECK_EQUAL(to_string(*cg_mol, Format::FASTA),
                          ">\nYTSLIEELIKKTEEQQKKNEEELKKLEEWAKKWNWF\n");
    }

    {
        // cyclic peptide
        std::unique_ptr<::RDKit::RWMol> mol(::RDKit::SmilesToMol(
            "CCC(C)[C@@H]1NC(=O)C2CCCN2C(=O)[C@H](CCCCN)NC(=O)[C@H](CC(N)=O)NC("
            "=O)[C@H](CC(C)C)NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](CC(C)C)"
            "NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H]2CCCN2C1=O"));
        bool use_residue_info = false;
        auto cg_mol = atomistic_to_cg(*mol, use_residue_info);
        BOOST_CHECK_EQUAL(to_string(*cg_mol, Format::HELM),
                          "PEPTIDE1{P.P.F.L.W.L.N.K.P.I}$PEPTIDE1,PEPTIDE1,10:"
                          "R2-1:R1$$$V2.0");
    }
}

BOOST_AUTO_TEST_CASE(Test_annotated_atomistic_to_cg)
{
    auto atomistic_mol =
        RDKit::v2::FileParsers::MolFromPDBFile(mmshare_testfile("1dng.pdb"));
    auto cg_mol = atomistic_to_cg(*atomistic_mol);

    BOOST_CHECK_EQUAL(to_string(*cg_mol, Format::HELM),
                      "PEPTIDE1{Q.A.P.A.Y.E.E.A.A.E.E.L.A.K.S}$$$$V2.0");
}

BOOST_AUTO_TEST_CASE(Test_cg_to_atomistic)
{
    // More tests validation tests are available in python, these are to ensure
    // memtest is ran
    {
        auto cg_mol =
            to_rdkit("PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0",
                     Format::HELM);
        auto atomistic_mol = cg_to_atomistic(*cg_mol);
        BOOST_CHECK_EQUAL(::RDKit::MolToSmiles(*atomistic_mol),
                          "C[C@@H]1NC(=O)[C@H](C)NC(=O)[C@@H](N)CSSC[C@@H](C(="
                          "O)O)NC(=O)[C@H](C)NC1=O");
    }
    {
        auto cg_mol = to_rdkit("PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,"
                               "PEPTIDE2,2:R3-1:R1$$$V2.0");
        auto atomistic_mol = cg_to_atomistic(*cg_mol);
        BOOST_CHECK_EQUAL(
            ::RDKit::MolToSmiles(*atomistic_mol),
            "N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)N[C@@H](CS)C(=O)N[C@@H](CCC(="
            "O)O)C(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)NCC(=O)O");
    }
    {
        // error condition -- invalid monomer
        auto cg_mol = to_rdkit("PEPTIDE1{F.Y.Z.G.R.L}$$$$V2.0");
        BOOST_CHECK_THROW(cg_to_atomistic(*cg_mol), std::out_of_range);
    }
    {
        // nucleic acids not yet supported
        auto cg_mol = to_rdkit("RNA1{A.C.G.T}$$$$V2.0");
        BOOST_CHECK_THROW(cg_to_atomistic(*cg_mol), std::runtime_error);
    }
}

BOOST_AUTO_TEST_CASE(Test_reordering_residues)
{
    ::RDKit::RWMol cg_mol;
    add_monomer(cg_mol, "D", 3, "PEPTIDE1");
    add_monomer(cg_mol, "C", 1, "PEPTIDE1");
    add_monomer(cg_mol, "B", 4, "PEPTIDE1");
    add_monomer(cg_mol, "A", 2, "PEPTIDE1");
    add_connection(cg_mol, 3, 2);
    add_connection(cg_mol, 2, 1);
    add_connection(cg_mol, 1, 0);
    assign_chains(cg_mol);

    BOOST_CHECK_EQUAL(to_string(cg_mol, Format::HELM),
                      "PEPTIDE1{A.B.C.D}$$$$V2.0");
}