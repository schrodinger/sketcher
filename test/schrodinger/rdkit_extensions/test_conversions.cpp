#define BOOST_TEST_MODULE test_rdkit_extensions_conversions

#include <fstream>
#include <iterator>
#include <string>
#include <tuple>
#include <vector>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "test_common.h"
#include "schrodinger/rdkit_extensions/atomistic_conversions.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/rdkit_extensions/helm.h"

using namespace schrodinger::rdkit_extensions;
namespace bdata = boost::unit_test::data;

using TestData = std::pair<std::string, std::string>;

namespace std
{
std::ostream& operator<<(std::ostream& os, const TestData& pair)
{
    os << pair.first << " " << pair.second;
    return os;
}
} // namespace std

// Example HELM strings that we should be able to roundtrip
// HELM -> atomistic -> monomeric -> HELM with the different residue
// identification methods (pdb info, substance groups, smarts matching)
static const std::vector<std::string> ROUNDTRIP_HELM_TEST_SET = {
    "PEPTIDE1{F.Y.K.A.R.L}$$$$V2.0",
    "PEPTIDE1{A.A.P.L}$$$$V2.0",
    "PEPTIDE1{A.K.A}$$$$V2.0",
    "PEPTIDE1{A.C.D.R.L}$$$$V2.0",
    // Examples with glycine, which requires a special query
    "PEPTIDE1{F.Y.K.G.R.L}$$$$V2.0",
    "PEPTIDE1{A.C.D.G.L}$$$$V2.0",
    // Cyclic peptide with R2-R1 closure
    "PEPTIDE1{A.A.A.A.A.A.A.A}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0",
    // Cyclic peptide using disulfide R3 - R3 bond
    "PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0",
    // Branches using R3 connections
    "PEPTIDE1{A.D(C)P}$$$$V2.0",
    "PEPTIDE1{A.C(A)G.C.D(K)A.D}$$$$V2.0",
    // Example with endcaps
    "PEPTIDE1{[ac].P.A.D(A)K.[am]}$$$$V2.0",
    // Example with a branch and a disulfide bond
    "PEPTIDE1{C.A.A.A.C.G.D(K)L.A}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0",
};

static boost::shared_ptr<RDKit::RWMol>
file_to_rdkit(const std::string& filename)
{
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open file: " + filename);
    }

    std::string content((std::istreambuf_iterator<char>(file)),
                        std::istreambuf_iterator<char>());
    file.close();
    return to_rdkit(content);
}

[[maybe_unused]] static std::unique_ptr<RDKit::ROMol>
resolve_his(const RDKit::ROMol& mol)
{
    // Some structures may contain different protonation states for histidine,
    // but we currently map all of them to the same single letter code 'H' in
    // the monomeric representation. Since we want to test the roundtrip
    // conversion against the original, we need to resolve the histidine
    // protonation state to what is in the monomer database.
    std::string smiles = RDKit::MolToSmiles(mol);

    std::vector<std::string> targets = {"Cc1c[nH]cn1"};
    std::string replace_with = "Cc1cnc[nH]1";
    for (const auto& target : targets) {
        size_t pos = 0;
        while ((pos = smiles.find(target, pos)) != std::string::npos) {
            smiles.replace(pos, target.length(), replace_with);
            pos += replace_with.length();
        }
    }
    return std::unique_ptr<RDKit::ROMol>(RDKit::SmilesToMol(smiles));
}

[[maybe_unused]] static void neutralize_molecule(RDKit::ROMol& mol)
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
}

[[maybe_unused]] static void remove_solvents(RDKit::RWMol& rwmol)
{
    std::vector<unsigned int> atoms_to_remove;
    for (const auto& atom : rwmol.atoms()) {
        const auto res_info = static_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());

        if (res_info != nullptr) {
            std::string residue_name = res_info->getResidueName();
            if (residue_name == "HOH" || residue_name == "S04") {
                atoms_to_remove.push_back(atom->getIdx());
            }
        }
    }
    rwmol.beginBatchEdit();
    for (unsigned int atom_idx : atoms_to_remove) {
        rwmol.removeAtom(atom_idx);
    }
    rwmol.commitBatchEdit();
}

BOOST_DATA_TEST_CASE(
    TestAtomisticSmilesToMonomeric,
    bdata::make(std::vector<TestData>{
        {"NCCCC[C@H](NC(=O)[C@H](CS)NC(=O)[C@@H](Cc1ccccc1)NC(=O)CNC(=O)[C@H]("
         "CCCNC(=N)N)NC(=O)[C@H](C)N)C(=O)N(C)[C@@H](C)C(=O)N[C@@H](CCC(=O)O)C("
         "=O)N[C@@H](CC(=O)O)C(=O)N[C@@H](C)C(=O)O",
         "PEPTIDE1{A.R.G.[dF].C.K.[meA].E.D.A}$$$$V2.0"},
        {"CC(C)C[C@@H]1NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H](CCCCN)NC(=O)["
         "C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H]2CCCN2C(=O)[C@H](C(C)O)NC("
         "=O)[C@H](CCC(=O)O)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](CC(N)=O)NC1=O",

         // The SMILEs is for a Threonine residue, which is missing one
         // chirality, so it doesn't match the template in the monomer DB.
         "PEPTIDE1{N.F.E.[CC(O)[C@H](N[*:1])C(=O)[*:2]].[dP]"
         ".P.F.K.W.L}$PEPTIDE1,PEPTIDE1,10:R2-1:R1$$$V2.0"}}),
    test_data)
{
    boost::shared_ptr<RDKit::ROMol> atomistic_mol(
        RDKit::SmilesToMol(test_data.first));

    bool try_residue_info = false;
    auto monomer_mol = toMonomeric(*atomistic_mol, try_residue_info);

    auto helm_result = to_string(*monomer_mol, Format::HELM);
    BOOST_TEST(helm_result == test_data.second);

    auto roundtrip_monomer_mol = to_rdkit(helm_result);
    auto roundtrip_atomistic = toAtomistic(*roundtrip_monomer_mol);

    RDKit::SubstructMatchParameters params;
    params.maxMatches = 1;
    auto match =
        RDKit::SubstructMatch(*roundtrip_atomistic, *atomistic_mol, params);
    BOOST_REQUIRE(!match.empty());
    BOOST_TEST(match[0].size() == atomistic_mol->getNumAtoms());
    BOOST_TEST((roundtrip_atomistic->getNumAtoms() -
                atomistic_mol->getNumAtoms()) <= 2);
}

BOOST_DATA_TEST_CASE(
    TestMonomericToAtomistic,
    bdata::make(std::vector<TestData>{
        {"PEPTIDE1{A.D(C)P}$$$$V2.0", "C[C@H](N)C(=O)N[C@@H](CC(=O)N[C@@H](CS)"
                                      "C(=O)O)C(=O)N1CCC[C@H]1C(=O)O"},
        {// Map numbers should work
         "PEPTIDE1{Q.R.F.[CC(C)(S)[C@@H](N[*:1])C(=O)[*:2]].T.G.H.F.G.G.L.Y.[O="
         "C([C@H]1CCCN1[*:1])[*:2]].[O=C([C@@H](CS)N[*:1])[*:2]].N.G.P}$$$$V2."
         "0",
         "CC(C)C[C@H](NC(=O)CNC(=O)CNC(=O)[C@H](Cc1ccccc1)NC(=O)[C@H](Cc1cnc["
         "nH]1)NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@H](Cc1ccccc1)NC(=O)[C@"
         "H](CCCNC(=N)N)NC(=O)[C@@H](N)CCC(N)=O)C(C)(C)S)[C@@H](C)O)C(=O)N[C@@"
         "H](Cc1ccc(O)cc1)C(=O)N1CCC[C@@H]1C(=O)N[C@H](CS)C(=O)N[C@@H](CC(N)=O)"
         "C(=O)NCC(=O)N1CCC[C@H]1C(=O)O"},
        {// RGroups should also be valid attachment point types, and R3-R3 bonds
         // should work
         "PEPTIDE1{[ac].C.F.F.[dW].K.T.F.[*N[C@H](C(=O)O)S* "
         "|$_R1;;;;;;;_R3$|]}$PEPTIDE1,PEPTIDE1,2:R3-9:R3$$$V2.0",
         "CC(=O)N[C@H]1CSS[C@@H](C(=O)O)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H]([C@@"
         "H](C)O)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H]("
         "Cc2ccccc2)NC(=O)[C@H](Cc2ccccc2)NC1=O"},
        {// Validated using pyPept (see LEGAL-1585)
         // R1-R2 connection -- this should technically be a single chain, but
         // still works
         "PEPTIDE1{A.A.A.A}|PEPTIDE2{P.P.P.P}$PEPTIDE1,PEPTIDE2,1:R1-4:R2$$$",
         "C[C@H](NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@H](C)NC(=O)[C@@H]1CCCN1C("
         "=O)[C@@H]1CCCN1C(=O)[C@@H]1CCCN1C(=O)[C@@H]1CCCN1)C(=O)O"},
        {// R3-R3 connection (disulfide bond)
         "PEPTIDE1{C.W.R.[dS].R.Y.[am]}|PEPTIDE2{C.W.R.[dS].R.Y.[am]}$PEPTIDE1,"
         "PEPTIDE2,1:R3-1:R3$$$",
         "N=C(N)NCCC[C@H](NC(=O)[C@H](Cc1c[nH]c2ccccc12)NC(=O)[C@@H](N)CSSC[C@"
         "H](N)C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@"
         "H](CO)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](Cc1ccc(O)cc1)C(N)=O)C(=O)"
         "N[C@H](CO)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](Cc1ccc(O)cc1)C(N)=O"},
        {// R3-R1 connection, from middle of PEPTIDE1
         "PEPTIDE1{D.E.F.G}|PEPTIDE2{C.E}$PEPTIDE1,PEPTIDE2,2:R3-1:R1$$$V2.0",
         "N[C@@H](CC(=O)O)C(=O)N[C@@H](CCC(=O)N[C@@H](CS)C(=O)N[C@@H](CCC(=O)O)"
         "C(=O)O)C(=O)N[C@@H](Cc1ccccc1)C(=O)NCC(=O)O"},
        {// Validated using pyPept (see LEGAL-1585), disulfide bond
         "PEPTIDE1{C.A.A.A.C}$PEPTIDE1,PEPTIDE1,1:R3-5:R3$$$V2.0",
         "C[C@@H]1NC(=O)[C@H](C)NC(=O)[C@@H](N)CSSC[C@@H](C(=O)O)NC(=O)[C@H](C)"
         "NC1=O"},
        {// Validated using pyPept, multiple disulfide bonds (slightly different
         // from pyPept due to different hydrogen placement on histidine)
         "PEPTIDE1{[dC].[dS].C.[dS].S.[dL].M.[dD].K.[dE].C.[dV].Y.[dF].[dC].H."
         "A.D.I.I.W}$PEPTIDE1,PEPTIDE1,3:R3-11:R3|PEPTIDE1,PEPTIDE1,1:R3-15:R3$"
         "$$",
         "CC[C@H](C)[C@H](NC(=O)[C@H](CC(=O)O)NC(=O)[C@H](C)NC(=O)[C@H](Cc1cnc["
         "nH]1)NC(=O)[C@H]1CSSC[C@@H](N)C(=O)N[C@H](CO)C(=O)N[C@H]2CSSC[C@H]("
         "NC(=O)[C@@H](CCC(=O)O)NC(=O)[C@H](CCCCN)NC(=O)[C@@H](CC(=O)O)NC(=O)["
         "C@H](CCSC)NC(=O)[C@@H](CC(C)C)NC(=O)[C@H](CO)NC(=O)[C@@H](CO)NC2=O)C("
         "=O)N[C@H](C(C)C)C(=O)N[C@@H](Cc2ccc(O)cc2)C(=O)N[C@H](Cc2ccccc2)C(=O)"
         "N1)C(=O)N[C@H](C(=O)N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O)[C@@H](C)CC"},
        {"PEPTIDE1{C.C.[dN].C.S.A.K.W.C.R.D.H.S.R.C.C.[am]}$PEPTIDE1,PEPTIDE1,"
         "1:R3-9:R3|PEPTIDE1,PEPTIDE1,2:R3-15:R3|PEPTIDE1,PEPTIDE1,4:R3-16:R3$$"
         "$",
         "C[C@@H]1NC(=O)[C@H](CO)NC(=O)[C@@H]2CSSC[C@@H](C(N)=O)NC(=O)[C@@H]"
         "3CSSC[C@H](NC(=O)[C@@H](N)CSSC[C@H](NC(=O)[C@H](Cc4c[nH]c5ccccc45)NC("
         "=O)[C@H](CCCCN)NC1=O)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@@H](CC(=O)O)C(="
         "O)N[C@@H](Cc1cnc[nH]1)C(=O)N[C@@H](CO)C(=O)N[C@@H](CCCNC(=N)N)C(=O)"
         "N3)C(=O)N[C@H](CC(N)=O)C(=O)N2"},
        {// Validated using pyPept, R2-R1 closure
         "PEPTIDE1{[dD].D.I.K.E.I.Y.D.P.G}$PEPTIDE1,PEPTIDE1,10:R2-1:R1$$$",
         "CC[C@H](C)[C@@H]1NC(=O)[C@H](CC(=O)O)NC(=O)[C@@H](CC(=O)O)NC(=O)CNC(="
         "O)[C@@H]2CCCN2C(=O)[C@H](CC(=O)O)NC(=O)[C@H](Cc2ccc(O)cc2)NC(=O)[C@H]"
         "([C@@H](C)CC)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H](CCCCN)NC1=O"},
        {"PEPTIDE1{A.A.G.F.P.V.F.F}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$",
         "CC(C)[C@@H]1NC(=O)[C@@H]2CCCN2C(=O)[C@H](Cc2ccccc2)NC(=O)CNC(=O)[C@H]"
         "(C)NC(=O)[C@H](C)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@H](Cc2ccccc2)NC1=O"},
        {// From HELM paper https://pubs.acs.org/doi/10.1021/ci3001925 + SMILES
         // canonicalized with RDKit
         "PEPTIDE1{A.R.C.A.A.K.T.C.D.A}$PEPTIDE1,PEPTIDE1,8:R3-3:R3$$$",
         "C[C@H](N)C(=O)N[C@@H](CCCNC(=N)N)C(=O)N[C@H]1CSSC[C@@H](C(=O)N[C@@H]("
         "CC(=O)O)C(=O)N[C@@H](C)C(=O)O)NC(=O)[C@H]([C@@H](C)O)NC(=O)[C@H]("
         "CCCCN)NC(=O)[C@H](C)NC(=O)[C@H](C)NC1=O"}}),
    test_data)
{
    auto monomer_mol = to_rdkit(test_data.first);
    auto atomistic_mol = toAtomistic(*monomer_mol);
    auto smiles_result = RDKit::MolToSmiles(*atomistic_mol);
    BOOST_TEST(smiles_result == test_data.second);
}

BOOST_DATA_TEST_CASE(TestRoundtripHelm, bdata::make(ROUNDTRIP_HELM_TEST_SET),
                     helm_str)
{
    // Test atomistic -> monomeric -> HELM using SMARTS matching method
    auto monomer_mol = to_rdkit(helm_str);
    auto atomistic_mol = toAtomistic(*monomer_mol);

    bool try_residue_info = false;
    auto mm_roundtrip = toMonomeric(*atomistic_mol, try_residue_info);

    auto helm_result = to_string(*mm_roundtrip, Format::HELM);
    BOOST_TEST(helm_result == helm_str);
}

BOOST_DATA_TEST_CASE(TestRoundtripHelmPDBResidueInfo,
                     bdata::make(ROUNDTRIP_HELM_TEST_SET), test_data)
{
    // Test atomistic -> monomeric -> HELM using existing PDB residue info
    auto monomer_mol = to_rdkit(test_data);
    auto atomistic_mol = toAtomistic(*monomer_mol);

    bool try_residue_info = true;
    auto mm_roundtrip = toMonomeric(*atomistic_mol, try_residue_info);

    auto helm_result = to_string(*mm_roundtrip, Format::HELM);
    BOOST_TEST(helm_result == test_data);
}

BOOST_DATA_TEST_CASE(TestRoundtripSupGroups,
                     bdata::make(ROUNDTRIP_HELM_TEST_SET), test_data)
{
    // Convert HELM to mol block and back to RDKit ROMol to ensure SUP groups
    // work as residue annotations
    auto monomer_mol = to_rdkit(test_data);
    auto atomistic_mol = toAtomistic(*monomer_mol);
    auto molblock = to_string(*atomistic_mol, Format::MDL_MOLV3000);
    boost::shared_ptr<RDKit::ROMol> mol_from_block(to_rdkit(molblock));

    bool try_residue_info = true;
    auto mm_from_sup_groups = toMonomeric(*mol_from_block, try_residue_info);
    auto helm_result = to_string(*mm_from_sup_groups, Format::HELM);
    BOOST_TEST(helm_result == test_data);
}

BOOST_AUTO_TEST_CASE(TestSetPdbInfo)
{
    // PDB residue information should be set correctly when converting from HELM
    // to atomistic
    std::string helm_str = "PEPTIDE1{A.A.A.A}|PEPTIDE2{P.P.P.P}$PEPTIDE1,"
                           "PEPTIDE2,1:R1-4:R2$$$V2.0";
    auto monomer_mol = to_rdkit(helm_str);
    auto atomistic_mol = toAtomistic(*monomer_mol);

    // Test first atom (should be from first chain, alanine)
    auto atom0 = atomistic_mol->getAtomWithIdx(0);
    auto info0 =
        static_cast<RDKit::AtomPDBResidueInfo*>(atom0->getMonomerInfo());
    BOOST_REQUIRE(info0 != nullptr);
    BOOST_TEST(info0->getResidueName() == "ALA");
    BOOST_TEST(info0->getResidueNumber() == 1);
    BOOST_TEST(info0->getChainId() == "A");
    BOOST_TEST(info0->getName() == " CB ");

    // Test atom from second chain (proline)
    auto atom34 = atomistic_mol->getAtomWithIdx(34);
    auto info34 =
        static_cast<RDKit::AtomPDBResidueInfo*>(atom34->getMonomerInfo());
    BOOST_REQUIRE(info34 != nullptr);
    BOOST_TEST(info34->getResidueName() == "PRO");
    BOOST_TEST(info34->getResidueNumber() == 2);
    BOOST_TEST(info34->getChainId() == "B");
    BOOST_TEST(info34->getName() == " N  ");
}

BOOST_AUTO_TEST_CASE(TestRnaHELMToAtomistic)
{
    std::string helm_str = "RNA1{P.[dR](A)P.[dR](C)}$$$$V2.0";
    auto monomer_mol = to_rdkit(helm_str);
    auto atomistic_mol = toAtomistic(*monomer_mol);

    auto smiles_result = RDKit::MolToSmiles(*atomistic_mol);
    std::string expected_smiles =
        "Nc1ccn([C@@H]2C[C@H](O)[C@@H](COP(=O)(O)O[C@H]3C[C@@H](n4cnc5c(N)"
        "ncnc54)O[C@@H]3COP(=O)(O)O)O2)c(=O)n1";
    BOOST_TEST(smiles_result == expected_smiles);
}

BOOST_AUTO_TEST_CASE(TestShared11152)
{
    auto atomistic_mol = file_to_rdkit(testfile_path("shared_11152.mae"));
    auto monomer_mol = toMonomeric(*atomistic_mol);

    auto helm_result = to_string(*monomer_mol, Format::HELM);
    std::string expected_helm =
        "PEPTIDE1{V.Q.L.Q.Q.S.G.G.E.L.A.K.P.G.A.S.V.K.V.S.C.K.A.S.G.Y.T.F.S.S."
        "F.W.M.H.W.V.R.Q.A.P.G.Q.G.L.E.W.I.G.Y.I.N.P.R.S.G.Y.T.E.Y.N.E.I.F.R.D."
        "K.A.T.M.T.T.D.T.S.T.S.T.A.Y.M.E.L.S.S.L.R.S.E.D.T.A.V.Y.Y.C.A.S.F.L.G."
        "R.G.A.M.D.Y.W.G.Q.G.T.T.V.T.V.S.S}|PEPTIDE2{E.I.Q.M.T.Q.S.P.S.S.L.S.A."
        "S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.I.S.N.Y.L.A.W.Y.Q.Q.K.P.G.K.A.P.K.L.L.I."
        "Y.Y.T.S.K.I.H.S.G.V.P.S.R.F.S.G.S.G.S.G.T.D.Y.T.F.T.I.S.T.N.I.S.L.Q.P."
        "E.D.I.A.T.Y.Y.C.Q.L.L.P.G.N.T.F.P.Y.T.F.G.V.A.N.G.T.K.V.A}$PEPTIDE2,"
        "PEPTIDE2,23:R3-91:R3|PEPTIDE1,PEPTIDE1,21:R3-95:R3$$$V2.0";
    BOOST_TEST(helm_result == expected_helm);
}

BOOST_DATA_TEST_CASE(
    TestCustomBonds,
    bdata::make(std::vector<TestData>{
        {"2n65.pdb",
         "PEPTIDE1{V.A.R.G.W.K.R.K.C.P.L.F.G.K.G.G}|PEPTIDE2{V.A.R.G.W.K.R.K.C."
         "P.L.F.G.K.G.G}$PEPTIDE1,PEPTIDE2,9:R3-9:R3$$$V2.0"},
        {"4qaf.pdb",
         "PEPTIDE1{E.E.I.Q.D.V.S.G.T.W.Y.L.K.A.M.T.V.D.V.G.A.L.R.C.L.A.G.S.V.I."
         "P.T.T.L.T.T.L.E.G.G.N.L.E.A.K.V.T.M.H.I.K.G.R.S.Q.E.V.K.A.V.L.S.K.T."
         "D.E.P.G.I.Y.T.A.I.G.G.I.H.V.A.K.I.G.R.S.H.V.K.D.H.Y.I.F.Y.S.E.G.C.L."
         "S.G.V.P.V.P.G.V.W.L.V.[CCCCCC[C@H]1C[C@H]1CCCCCCCCCC(O)O.[*:1]]}|"
         "PEPTIDE2{V.S.G.T.W.Y.L.K.A.M.T.V.D.V.G.A.L.R.C.L.A.G.S.V.I.P.T.T.L.T."
         "T.N.L.E.A.K.V.T.E.V.K.A.V.L.S.K.T.D.E.P.G.I.Y.T.A.I.G.G.I.H.V.A.K.I."
         "G.R.S.D.H.Y.I.F.Y.S.E.G.C.L.S.G.V.P.V.P.G.V.W.L.V.[CCCCCC[C@H]1C[C@H]"
         "1CCCCCCCCCC(O)O.[*:1]]}|PEPTIDE3{E.V.V.K.F.M.D.V.Y.Q.R.S.Y.C.H.P.I.E."
         "T.L.V.D.I.F.Q.E.Y.P.D.E.I.E.Y.I.F.K.P.S.C.V.P.L.M.R.C.G.G.C.C.N.D.E."
         "G.L.E.C.V.P.T.E.E.S.N.I.T.M.Q.I.M.R.I.K.P.H.Q.G.Q.H.I.G.E.M.S.F.L.Q."
         "H.N.K.C.E.C.R.P.K.K.[OS(O)(O)O.[*:1].[*:2]].[OS(O)(O)O.[*:1]]}|"
         "PEPTIDE4{H.E.V.V.K.F.M.D.V.Y.Q.R.S.Y.C.H.P.I.E.T.L.V.D.I.F.Q.E.Y.P.D."
         "E.I.E.Y.I.F.K.P.S.C.V.P.L.M.R.C.G.G.C.C.N.D.E.G.L.E.C.V.P.T.E.E.S.N."
         "I.T.M.Q.I.M.R.I.K.P.H.Q.G.Q.H.I.G.E.M.S.F.L.Q.H.N.K.C.E.C.R.P.[OS(O)("
         "O)O.[*:1].[*:2]].[CC(O)O.[*:1]]}$PEPTIDE1,PEPTIDE1,24:R3-97:R3|"
         "PEPTIDE3,PEPTIDE4,39:R3-49:R3|PEPTIDE3,PEPTIDE3,45:R3-90:R3|PEPTIDE3,"
         "PEPTIDE4,48:R3-40:R3|PEPTIDE3,PEPTIDE3,49:R3-92:R3|PEPTIDE4,PEPTIDE4,"
         "46:R3-91:R3|PEPTIDE4,PEPTIDE4,50:R3-93:R3$$$V2.0"},
        {"5vav.pdb", "PEPTIDE1{G.R.C.T.Q.A.W.P.P.I.C.F.P.D}$PEPTIDE1,PEPTIDE1,"
                     "3:R3-11:R3|PEPTIDE1,PEPTIDE1,14:R2-1:R1$$$V2.0"},
        {"shared_11152.mae",
         "PEPTIDE1{V.Q.L.Q.Q.S.G.G.E.L.A.K.P.G.A.S.V.K.V.S.C.K.A.S.G.Y.T.F.S.S."
         "F.W.M.H.W.V.R.Q.A.P.G.Q.G.L.E.W.I.G.Y.I.N.P.R.S.G.Y.T.E.Y.N.E.I.F.R."
         "D.K.A.T.M.T.T.D.T.S.T.S.T.A.Y.M.E.L.S.S.L.R.S.E.D.T.A.V.Y.Y.C.A.S.F."
         "L.G.R.G.A.M.D.Y.W.G.Q.G.T.T.V.T.V.S.S}|PEPTIDE2{E.I.Q.M.T.Q.S.P.S.S."
         "L.S.A.S.V.G.D.R.V.T.I.T.C.R.A.S.Q.D.I.S.N.Y.L.A.W.Y.Q.Q.K.P.G.K.A.P."
         "K.L.L.I.Y.Y.T.S.K.I.H.S.G.V.P.S.R.F.S.G.S.G.S.G.T.D.Y.T.F.T.I.S.T.N."
         "I.S.L.Q.P.E.D.I.A.T.Y.Y.C.Q.L.L.P.G.N.T.F.P.Y.T.F.G.V.A.N.G.T.K.V.A}$"
         "PEPTIDE2,PEPTIDE2,23:R3-91:R3|PEPTIDE1,PEPTIDE1,21:R3-95:R3$$$V2.0"},
        {"1dng.pdb", "PEPTIDE1{Q.A.P.A.Y.E.E.A.A.E.E.L.A.K.S}$$$$V2.0"},
        {"ionized_residue.pdb", "PEPTIDE1{A.[OC(O)C[C@H](N[*:1])C(O)[*:2]].A.["
                                "OC(O)CC[C@H](N[*:1])C(O)[*:2]].A}$$$$V2.0"},
        {"AF-Q6Q6F6-F1-model_v1.pdb",
         "PEPTIDE1{M.S.S.R.E.R.R.S.D.L.Y.I.K.A.E.P.S.S.P.E.G.G.G.G.G.G.G.G.R.T."
         "S.P.G.G.A.S.S.D.S.S.Q.S.G.G.G.G.S.R.G.E.G.A.G.R.Y.S.P.P.L.Y.T.P.A.L."
         "R.C.H.F.K.E.E.G.A.D.G.A.E.E.G.S.T.G.S.G.G.G.R.C.K.Y.A.L.S.T.L.P.K.R."
         "L.C.L.V.C.G.D.V.A.S.G.Y.H.Y.G.V.A.S.C.E.A.C.K.A.F.F.K.R.T.I.Q.G.N.I."
         "E.Y.S.C.P.A.S.N.E.C.E.I.T.K.R.R.R.K.A.C.Q.A.C.R.F.T.K.C.L.K.V.G.M.L."
         "K.E.G.V.R.L.D.R.V.R.G.G.R.Q.K.Y.K.R.R.P.E.V.E.N.A.T.Y.Q.S.A.P.I.P.L."
         "R.K.E.G.E.K.G.S.S.S.I.I.V.S.H.L.L.V.A.E.P.E.K.L.F.A.M.P.D.P.L.Q.P.D."
         "T.A.Q.R.T.L.T.T.L.C.D.L.A.D.R.E.L.V.V.I.I.G.W.A.K.H.I.P.G.F.L.S.L.S."
         "L.A.D.Q.M.S.V.L.Q.S.V.W.L.E.V.L.V.L.G.V.A.Y.R.S.L.G.C.E.D.E.V.V.F.A."
         "E.D.F.V.L.D.E.E.M.S.R.V.A.G.L.T.E.L.N.A.A.I.S.Q.L.A.R.R.F.R.A.L.Q.L."
         "D.R.E.E.F.V.M.L.K.A.I.A.L.T.N.S.D.S.V.Y.I.E.D.M.E.A.V.Q.K.L.R.D.L.L."
         "H.Q.A.L.L.E.L.E.V.Q.R.R.P.D.D.P.Q.R.A.G.R.L.L.L.T.L.P.L.L.R.Q.T.A.G."
         "R.A.L.T.T.F.Y.S.I.K.T.R.G.G.V.P.M.H.K.L.F.L.E.M.L.E.A.M.M.D.S.P}$$$$"
         "V2.0"}}),
    test_data)
{
    // Test that atomistic -> monomeric -> atomistic produces the same structure
    // as the original atomistic structure
    auto file_path = testfile_path(test_data.first);
    auto original_atomistic = file_to_rdkit(file_path);
    remove_solvents(*original_atomistic);

    auto monomer_mol = toMonomeric(*original_atomistic);
    auto helm_result = to_string(*monomer_mol, Format::HELM);
    BOOST_TEST(helm_result == test_data.second);

#ifndef WIN32
    // Substructure matching in this test causes segv on windows
    auto roundtrip_atomistic = toAtomistic(*monomer_mol);

    // Process these to ensure substructure matching works to verify equivalence
    auto original_atomistic_his = resolve_his(*original_atomistic);
    auto roundtrip_atomistic_his = resolve_his(*roundtrip_atomistic);
    neutralize_molecule(*roundtrip_atomistic_his);
    neutralize_molecule(*original_atomistic_his);

    // Verify the roundtrip produces a molecule with expected number of atoms,
    // roundtripped structure has additional oxygen atoms at the end of the
    // chains so there may be some extra atoms in the roundtrip (up to 2 *
    // number of chains)
    RDKit::SubstructMatchParameters params;
    params.maxMatches = 1;
    auto match = RDKit::SubstructMatch(*roundtrip_atomistic_his,
                                       *original_atomistic_his, params);
    BOOST_REQUIRE(!match.empty());
    BOOST_TEST(match[0].size() == original_atomistic_his->getNumAtoms());
    BOOST_TEST((roundtrip_atomistic_his->getNumAtoms() -
                original_atomistic_his->getNumAtoms()) <= 8);
#endif
}

BOOST_AUTO_TEST_CASE(TestInconsistentCycleStartPoint)
{
    auto roundtrip_helm = [](const std::string& helm_str) {
        auto monomer_mol = to_rdkit(helm_str);
        auto atomistic_mol = toAtomistic(*monomer_mol);
        bool try_residue_info = false;
        auto mm_roundtrip = toMonomeric(*atomistic_mol, try_residue_info);
        return to_string(*mm_roundtrip, Format::HELM);
    };

    // Atomistic -> Monomeric conversion can lead to inconsistent cycle start
    // points since there is no obvious place to start the cycle or a
    // canonicalization method. See SHARED-11738
    auto helm1 = "PEPTIDE1{A.F.P.V.F.F.A.A}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0";
    auto helm2 = "PEPTIDE1{F.P.V.F.F.A.A.A}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0";
    auto helm3 = "PEPTIDE1{P.V.F.F.A.A.A.F}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0";
    BOOST_TEST(roundtrip_helm(helm1) == helm2);
    BOOST_TEST(roundtrip_helm(helm2) == helm3);
    BOOST_TEST(roundtrip_helm(helm3) ==
               "PEPTIDE1{V.F.F.A.A.A.F.P}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestAutoDetectIfNoResidueInfo)
{
    {
        // If there is no residue information, the atomistic -> monomeric
        // conversion will fall back to the SMARTS matching method
        auto smiles = "CC(C)C[C@@H]1NC(=O)[C@H](Cc2c[nH]c3ccccc23)NC(=O)[C@H]("
                      "CCCCN)NC(=O)[C@H](Cc2ccccc2)NC(=O)[C@@H]2CCCN2C(=O)[C@H]"
                      "2CCCN2C(=O)[C@H](C(C)O)NC(=O)[C@H](CCC(=O)O)NC(=O)[C@H]("
                      "Cc2ccccc2)NC(=O)[C@H](CC(N)=O)NC1=O";
        boost::shared_ptr<RDKit::ROMol> atomistic_mol(
            RDKit::SmilesToMol(smiles));
        auto monomer_mol = toMonomeric(*atomistic_mol);
        auto helm_result = to_string(*monomer_mol, Format::HELM);

        // The SMILEs is for a Threonine residue, which is missing one
        // chirality, so it doesn't match the template in the monomer DB.
        BOOST_TEST(helm_result ==
                   "PEPTIDE1{N.F.E.[CC(O)[C@H](N[*:1])C(=O)[*:2]]"
                   ".[dP].P.F.K.W.L}$PEPTIDE1,PEPTIDE1,10:R2-1:R1$$$V2.0");
    }

    {
        // If we have a structure that doesn't resemble a sequence of amino
        // acids at all, it will be placed into a single SMILES monomer
        auto smiles = "CC(C)c1ccccc1";
        boost::shared_ptr<RDKit::ROMol> atomistic_mol(
            RDKit::SmilesToMol(smiles));
        auto monomer_mol = toMonomeric(*atomistic_mol);
        auto helm_result = to_string(*monomer_mol, Format::HELM);
        BOOST_TEST(helm_result == "CHEM1{[CC(C)c1ccccc1]}$$$$V2.0");
    }

    {
        // But if we have a single known amino acid, it will be identified as
        // such
        auto smiles = "N[C@@H](CC(=O)O)C(=O)O";
        boost::shared_ptr<RDKit::ROMol> atomistic_mol(
            RDKit::SmilesToMol(smiles));
        auto monomer_mol = toMonomeric(*atomistic_mol);
        auto helm_result = to_string(*monomer_mol, Format::HELM);
        BOOST_TEST(helm_result == "PEPTIDE1{D}$$$$V2.0");
    }
}

BOOST_AUTO_TEST_CASE(TestNoAtomsToMonomeric)
{
    // Test that we handle empty molecules correctly
    boost::shared_ptr<RDKit::ROMol> atomistic_mol(new RDKit::ROMol());
    BOOST_CHECK_THROW(toMonomeric(*atomistic_mol), std::runtime_error);
}

BOOST_DATA_TEST_CASE(
    TestAtomisticToMonomericSmilesMonomers,
    bdata::make(std::vector<std::string>{
        // Use of disulfide bond in SMILES monomer, matches standard backbone
        "PEPTIDE1{[ac].[CC(C)(S[*:3])[C@@H](N[*:1])C(=O)[*:2]].I.P.R.G.D.["
        "COc1ccc(C[C@H](N[*:1])C(=O)[*:2])cc1].R.C.[am]}$PEPTIDE1,PEPTIDE1,2:"
        "R3-10:R3$$$V2.0",
        // ----
        "PEPTIDE1{D.[CC(C)(S[*:3])[C@@H](N[*:1])C(=O)[*:2]].F.W.[NCCC[C@@H](N[*"
        ":1])C(=O)[*:2]].Y.C.V}$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$$V2.0",
        // ----
        "PEPTIDE1{F.W.[CC(C)(C[*:3])[C@@H](N[*:1])C(=O)[*:2]].P.A.G.C.K}$"
        "PEPTIDE1,PEPTIDE1,3:R3-7:R3$$$V2.0",
        // Example of R3 attachment point use on non-standard backbone (doesn't
        // match general AA query)
        "PEPTIDE1{F.W.[CC(C)(C[*:3])C(CO[*:1])C(=N)[*:2]].P.A.G.C.K}$PEPTIDE1,"
        "PEPTIDE1,3:R3-7:R3$$$V2.0",
        // Disulfide bond on a SMILES monomer not identified by the cysteine
        // query
        "PEPTIDE1{[ac].[CC(C)(CS[*:3])[C@@H](N[*:1])C(=O)[*:2]].I.P.R.G.D.["
        "COc1ccc(C[C@H](N[*:1])C(=O)[*:2])cc1].R.C.[am]}$PEPTIDE1,PEPTIDE1,2:"
        "R3-10:R3$$$V2.0",
        // Various map number uses, but ensure that the sulfur is not given an
        // R3 attachment point since it isn't used
        "PEPTIDE1{Q.R.F.[CC(C)(S)[C@@H](N[*:1])C(=O)[*:2]].T.G.H.F.G.G.L.Y.[O="
        "C([C@H]1CCCCCN1[*:1])[*:2]].[O=C([C@@H](CCS)N[*:1])[*:2]].N.G.P}$$$$"
        "V2.0"}),
    helm_str)
{
    // Test SMILES -> MonomerMol where some of the monomers are SMILES monomers
    // and the correct attachment points should be identified and used in their
    // connections. Test by roundtripping through HELM -> monomeric -> atomistic
    // -> monomeric -> ROUNDTRIP_HELM; HELM and ROUNDTRIP_HELM should be
    // identical
    auto monomer_mol = to_rdkit(helm_str);
    auto atomistic_mol = toAtomistic(*monomer_mol);
    bool use_residue_info = false;
    auto roundtrip_monomer_mol = toMonomeric(*atomistic_mol, use_residue_info);
    auto helm_sountrip = to_string(*roundtrip_monomer_mol, Format::HELM);
    BOOST_TEST(helm_sountrip == helm_str);
}

BOOST_AUTO_TEST_CASE(TestShared11819)
{
    auto monomer_mol = to_rdkit("CHEM1{[c1ccccc1]}$$$$V2.0");

    bool is_smiles = false;
    BOOST_TEST(monomer_mol->getAtomWithIdx(0)->getPropIfPresent(SMILES_MONOMER,
                                                                is_smiles));
    BOOST_TEST(is_smiles);

    auto atomistic_mol = toAtomistic(*monomer_mol);
    BOOST_TEST(RDKit::MolToSmiles(*atomistic_mol) == "c1ccccc1");
}
