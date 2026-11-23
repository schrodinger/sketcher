#define BOOST_TEST_MODULE test_rdkit_extensions_conversions

#include <algorithm>
#include <fstream>
#include <cmath>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <numbers>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/shared_ptr.hpp>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Depictor/DepictUtils.h>
#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "test_common.h"
#include "schrodinger/rdkit_extensions/atomistic_conversions.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

using namespace schrodinger::rdkit_extensions;
namespace bdata = boost::unit_test::data;

using TestData = std::pair<std::string, std::string>;

struct ResetMonomerDatabase {
    ~ResetMonomerDatabase()
    {
        MonomerDatabase::instance().resetMonomerDefinitions();
    }
};

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
    "PEPTIDE1{C}|PEPTIDE2{A.D.P}$PEPTIDE2,PEPTIDE1,2:R3-1:R1$$$V2.0",
    "PEPTIDE1{A.C.G.C.D.A.D}|PEPTIDE2{A}|PEPTIDE3{K}$PEPTIDE1,PEPTIDE2,2:R3-1:"
    "R1|PEPTIDE1,PEPTIDE3,5:R3-1:R1$$$V2.0",
    // Example with endcaps
    "PEPTIDE1{[ac].P.A.D.K.[am]}|PEPTIDE2{A}$PEPTIDE1,PEPTIDE2,4:R3-1:R1$$$V2."
    "0",
    // Example with a branch and a disulfide bond
    "PEPTIDE1{C.A.A.A.C.G.D.L.A}|PEPTIDE2{K}$PEPTIDE1,PEPTIDE1,1:R3-5:R3|"
    "PEPTIDE1,PEPTIDE2,7:R3-1:R1$$$V2.0",
    // Disulfide bond between two chains
    "PEPTIDE1{C}|PEPTIDE2{C}$PEPTIDE1,PEPTIDE2,1:R3-1:R3$$$V2.0",
};

BOOST_AUTO_TEST_CASE(TestAtomisticCoordinatesUsingMonomerLayout)
{
    // A-G-S checks that atomistic layout preserves the graph, follows the
    // monomer order, and still accepts the explicit coord-map API.
    auto monomer_mol = helm::helm_to_rdkit("PEPTIDE1{A.G.S}$$$$V2.0");
    auto atomistic_mol = toAtomistic(*monomer_mol);
    const auto num_atoms = atomistic_mol->getNumAtoms();
    const auto num_bonds = atomistic_mol->getNumBonds();

    compute2DCoords(*atomistic_mol);
    BOOST_TEST(atomistic_mol->getNumConformers() == 1);
    BOOST_TEST(atomistic_mol->getNumAtoms() == num_atoms);
    BOOST_TEST(atomistic_mol->getNumBonds() == num_bonds);

    std::map<int, std::vector<RDGeom::Point3D>> residue_positions;
    const auto& conf = atomistic_mol->getConformer();
    for (const auto* atom : atomistic_mol->atoms()) {
        const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        BOOST_REQUIRE(info != nullptr);
        const auto p = conf.getAtomPos(atom->getIdx());
        BOOST_TEST(std::isfinite(p.x));
        BOOST_TEST(std::isfinite(p.y));
        residue_positions[info->getResidueNumber()].push_back(p);
    }
    BOOST_REQUIRE(residue_positions.size() == 3);
    // The monomer layout is a left-to-right chain. The atomistic residues
    // must follow that order, even though zigzag bounds can overlap in x.
    double previous_center_x = -std::numeric_limits<double>::infinity();
    for (const auto& [number, positions] : residue_positions) {
        double min_x = std::numeric_limits<double>::infinity();
        double max_x = -min_x;
        for (const auto& p : positions) {
            min_x = std::min(min_x, p.x);
            max_x = std::max(max_x, p.x);
        }
        const auto center_x = (min_x + max_x) / 2;
        BOOST_TEST(center_x > previous_center_x);
        previous_center_x = center_x;
    }

    BOOST_CHECK_NO_THROW(compute2DCoords(*atomistic_mol, {0}));
}

BOOST_AUTO_TEST_CASE(TestPeptideBackboneZigzag)
{
    // Four alanines isolate the straight zigzag: regular backbone bonds,
    // shorter C=O and side-chain bonds, and alternating side-chain sides.
    auto monomer_mol = helm::helm_to_rdkit("PEPTIDE1{A.A.A.A}$$$$V2.0");
    auto atomistic_mol = toAtomistic(*monomer_mol);
    compute2DCoords(*atomistic_mol);

    struct BackbonePositions {
        RDGeom::Point3D nitrogen, alpha, carbonyl, oxygen, beta;
    };
    std::map<int, BackbonePositions> positions;
    const auto& conf = atomistic_mol->getConformer();
    for (const auto* atom : atomistic_mol->atoms()) {
        const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        BOOST_REQUIRE(info != nullptr);
        auto name = info->getName();
        name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
        const auto p = conf.getAtomPos(atom->getIdx());
        auto& backbone = positions[info->getResidueNumber()];
        if (name == "N")
            backbone.nitrogen = p;
        else if (name == "CA")
            backbone.alpha = p;
        else if (name == "C")
            backbone.carbonyl = p;
        else if (name == "O")
            backbone.oxygen = p;
        else if (name == "CB")
            backbone.beta = p;
    }
    BOOST_REQUIRE(positions.size() == 4);
    const auto bond_length = RDDepict::BOND_LEN;
    const auto sidechain_bond_length = RDDepict::BOND_LEN * 0.7;
    double previous_sidechain_dy = 0;
    for (int residue_number = 1; residue_number <= 4; ++residue_number) {
        const auto& current = positions.at(residue_number);
        const auto n_ca = current.alpha - current.nitrogen;
        const auto ca_c = current.carbonyl - current.alpha;
        const auto c_o = current.oxygen - current.carbonyl;
        const auto ca_cb = current.beta - current.alpha;
        for (const auto& bond : {n_ca, ca_c})
            BOOST_TEST(std::abs(std::hypot(bond.x, bond.y) - bond_length) <
                       1e-5);
        BOOST_TEST(std::abs(std::hypot(c_o.x, c_o.y) - sidechain_bond_length) <
                   1e-5);
        BOOST_TEST(std::abs(std::hypot(ca_cb.x, ca_cb.y) -
                            sidechain_bond_length) < 1e-5);
        BOOST_TEST(n_ca.x > 0);
        BOOST_TEST(ca_c.x > 0);
        BOOST_TEST(std::abs(std::atan2(std::abs(n_ca.y), n_ca.x) -
                            std::numbers::pi / 6) < 1e-5);
        BOOST_TEST(std::abs(std::atan2(std::abs(ca_c.y), ca_c.x) -
                            std::numbers::pi / 6) < 1e-5);
        BOOST_TEST(n_ca.y * ca_c.y < 0);
        const auto ca_n = current.nitrogen - current.alpha;
        const auto cosine_at_alpha =
            (ca_n.x * ca_c.x + ca_n.y * ca_c.y) / (bond_length * bond_length);
        BOOST_TEST(std::abs(cosine_at_alpha + 0.5) < 1e-5);
        BOOST_TEST(std::abs(c_o.x) < 1e-5);
        BOOST_TEST(c_o.y * n_ca.y < 0);
        BOOST_TEST(std::abs(ca_cb.x) < 1e-5);
        BOOST_TEST(ca_cb.y * n_ca.y > 0);
        if (residue_number > 1)
            BOOST_TEST(ca_cb.y * previous_sidechain_dy < 0);
        previous_sidechain_dy = ca_cb.y;
        if (residue_number < 4) {
            const auto peptide =
                positions.at(residue_number + 1).nitrogen - current.carbonyl;
            BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) -
                                bond_length) < 1e-5);
            BOOST_TEST(peptide.x > 0);
            BOOST_TEST(std::abs(std::atan2(std::abs(peptide.y), peptide.x) -
                                std::numbers::pi / 6) < 1e-5);
            BOOST_TEST(peptide.y * n_ca.y > 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestStraightChainProlineGeometry)
{
    // Check P and dP separately: N-CA opens horizontally while the incoming
    // peptide link stays trans and the ring faces the preceding side chain.
    for (const auto* helm :
         {"PEPTIDE1{A.P.A}$$$$V2.0", "PEPTIDE1{A.[dP].A}$$$$V2.0"}) {
        auto monomer_mol = helm::helm_to_rdkit(helm);
        auto mol = toAtomistic(*monomer_mol);
        compute2DCoords(*mol);
        std::map<int, std::map<std::string, RDGeom::Point3D>> positions;
        for (const auto* atom : mol->atoms()) {
            const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                atom->getMonomerInfo());
            if (!info)
                continue;
            auto name = info->getName();
            name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
            positions[info->getResidueNumber()][name] =
                mol->getConformer().getAtomPos(atom->getIdx());
        }
        const auto& before = positions.at(1);
        const auto& proline = positions.at(2);
        BOOST_TEST_CONTEXT(helm)
        {
            BOOST_TEST(std::abs(proline.at("N").y - proline.at("CA").y) < 1e-5);
            BOOST_TEST(proline.at("CA").x > proline.at("N").x);
            BOOST_TEST((before.at("CB").y - before.at("CA").y) *
                           (proline.at("CB").y - proline.at("CA").y) >
                       0);
            const auto peptide_bond = proline.at("N") - before.at("C");
            const auto before_alpha = before.at("CA") - before.at("C");
            const auto proline_alpha = proline.at("CA") - proline.at("N");
            const auto cross = [&](const RDGeom::Point3D& p) {
                return peptide_bond.x * p.y - peptide_bond.y * p.x;
            };
            BOOST_TEST(cross(before_alpha) * cross(proline_alpha) < 0);
        }
    }
}

BOOST_AUTO_TEST_CASE(TestDefaultCoordinateLayoutWithoutResidueInfo)
{
    // A small molecule without residue metadata must use ordinary 2D layout.
    std::unique_ptr<RDKit::RWMol> mol(RDKit::SmilesToMol("CCO"));
    compute2DCoords(*mol);
    BOOST_TEST(mol->getNumConformers() == 1);
}

BOOST_AUTO_TEST_CASE(TestUnsupportedPolymersUseRdkitCoordinates)
{
    // Cyclic peptide backbones have no starting residue for the guided walk;
    // RNA has no peptide backbone. Both use RDKit's complete graph layout.
    for (const auto& helm :
         {"PEPTIDE1{A.A.A.A.A.A.A.A}$PEPTIDE1,PEPTIDE1,8:R2-1:R1$$$V2.0",
          "RNA1{P.[dR](A)P.[dR](C)}$$$$V2.0"}) {
        auto monomers = helm::helm_to_rdkit(helm);
        auto actual = toAtomistic(*monomers);
        auto expected = toAtomistic(*monomers);
        compute2DCoords(*actual);
        compute2DCoords(*expected, {}, false);
        const auto& actual_conf = actual->getConformer();
        const auto& expected_conf = expected->getConformer();
        BOOST_REQUIRE(actual->getNumAtoms() == expected->getNumAtoms());
        for (unsigned int i = 0; i < actual->getNumAtoms(); ++i) {
            const auto delta =
                actual_conf.getAtomPos(i) - expected_conf.getAtomPos(i);
            BOOST_TEST_CONTEXT(helm)
            {
                BOOST_TEST(std::hypot(delta.x, delta.y) < 1e-5);
            }
        }
    }
}

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

BOOST_AUTO_TEST_CASE(TestMonomerLayoutFromPdb)
{
    // 1dng has explicit H, proline, straight arms, and a turn. Check their
    // bond scales and orientations after atomistic placement from PDB data.
    std::ifstream file(testfile_path("1dng.pdb"));
    BOOST_REQUIRE(file.is_open());
    const std::string pdb((std::istreambuf_iterator<char>(file)),
                          std::istreambuf_iterator<char>());
    auto mol = to_rdkit(pdb, Format::PDB);
    const auto num_atoms = mol->getNumAtoms();
    compute2DCoords(*mol);
    BOOST_TEST(mol->getNumAtoms() == num_atoms);
    BOOST_TEST(mol->getNumConformers() == 1);
    for (unsigned int i = 0; i < num_atoms; ++i) {
        const auto p = mol->getConformer().getAtomPos(i);
        BOOST_TEST(std::isfinite(p.x));
        BOOST_TEST(std::isfinite(p.y));
    }
    std::map<int, std::map<std::string, RDGeom::Point3D>> named_positions;
    for (const auto* atom : mol->atoms()) {
        const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        if (!info)
            continue;
        auto name = info->getName();
        name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
        named_positions[info->getResidueNumber()][name] =
            mol->getConformer().getAtomPos(atom->getIdx());
    }
    const auto bond_length = RDDepict::BOND_LEN;
    const auto sidechain_bond_length = RDDepict::BOND_LEN * 0.7;
    unsigned int horizontal_residues = 0;
    unsigned int checked_amide_hydrogens = 0;
    for (const auto& [number, atoms] : named_positions) {
        if (!atoms.contains("N") || !atoms.contains("CA") ||
            !atoms.contains("C") || !atoms.contains("O"))
            continue;
        const auto n_ca = atoms.at("CA") - atoms.at("N");
        const auto ca_c = atoms.at("C") - atoms.at("CA");
        const auto c_o = atoms.at("O") - atoms.at("C");
        for (const auto& bond : {n_ca, ca_c})
            BOOST_TEST(std::abs(std::hypot(bond.x, bond.y) - bond_length) <
                       1e-5);
        BOOST_TEST(std::abs(std::hypot(c_o.x, c_o.y) - sidechain_bond_length) <
                   1e-5);
        if (atoms.contains("OXT")) {
            const auto c_oxt = atoms.at("OXT") - atoms.at("C");
            BOOST_TEST(std::abs(std::hypot(c_oxt.x, c_oxt.y) -
                                sidechain_bond_length) < 1e-5);
        }
        if (atoms.contains("CB")) {
            const auto ca_cb = atoms.at("CB") - atoms.at("CA");
            BOOST_TEST(std::abs(std::hypot(ca_cb.x, ca_cb.y) -
                                (number == 3 ? bond_length
                                             : sidechain_bond_length)) < 1e-5);
        }
        const auto backbone = atoms.at("C") - atoms.at("N");
        if (atoms.contains("HA") && atoms.contains("CB")) {
            const auto ca_ha = atoms.at("HA") - atoms.at("CA");
            const auto ca_cb = atoms.at("CB") - atoms.at("CA");
            BOOST_TEST(std::abs(std::hypot(ca_ha.x, ca_ha.y) -
                                sidechain_bond_length) < 1e-5);
            BOOST_TEST(std::abs(ca_ha.x * backbone.x + ca_ha.y * backbone.y) <
                       1e-5);
            BOOST_TEST(ca_ha.x * ca_cb.x + ca_ha.y * ca_cb.y < 0);
        }
        const auto previous = named_positions.find(number - 1);
        if (atoms.contains("H") && previous != named_positions.end() &&
            previous->second.contains("C")) {
            const auto n_h = atoms.at("H") - atoms.at("N");
            const auto n_c = previous->second.at("C") - atoms.at("N");
            const auto n_ca_length = std::hypot(n_ca.x, n_ca.y);
            const auto n_c_length = std::hypot(n_c.x, n_c.y);
            const auto bisector = n_ca / n_ca_length + n_c / n_c_length;
            const auto bisector_length = std::hypot(bisector.x, bisector.y);
            const auto expected =
                bisector * (-sidechain_bond_length / bisector_length);
            BOOST_TEST(std::hypot(n_h.x - expected.x, n_h.y - expected.y) <
                       1e-5);
            ++checked_amide_hydrogens;
        }
        if (std::abs(backbone.y) < 1e-5) {
            ++horizontal_residues;
            BOOST_TEST(std::abs(c_o.x) < 1e-5);
            BOOST_TEST(std::abs(std::atan2(std::abs(n_ca.y), std::abs(n_ca.x)) -
                                std::numbers::pi / 6) < 1e-5);
            BOOST_TEST(std::abs(std::atan2(std::abs(ca_c.y), std::abs(ca_c.x)) -
                                std::numbers::pi / 6) < 1e-5);
        }
        const auto following = named_positions.find(number + 1);
        if (following != named_positions.end() &&
            following->second.contains("N")) {
            const auto peptide = following->second.at("N") - atoms.at("C");
            BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) -
                                bond_length) < 1e-5);
        }
    }
    BOOST_TEST(horizontal_residues >= 8);
    BOOST_TEST(checked_amide_hydrogens >= 8);
    BOOST_TEST(named_positions.at(3).at("N").x >
               named_positions.at(1).at("N").x);
    BOOST_TEST(named_positions.at(13).at("N").x >
               named_positions.at(15).at("N").x);
    const auto peptide_side = [&](int previous, int next) {
        const auto& before = named_positions.at(previous);
        const auto& after = named_positions.at(next);
        const auto bond = after.at("N") - before.at("C");
        const auto before_alpha = before.at("CA") - before.at("C");
        const auto after_alpha = after.at("CA") - after.at("N");
        const auto cross = [&](const RDGeom::Point3D& vector) {
            return bond.x * vector.y - bond.y * vector.x;
        };
        return cross(before_alpha) * cross(after_alpha);
    };
    BOOST_TEST(peptide_side(1, 2) < 0); // trans
    BOOST_TEST(peptide_side(2, 3) < 0); // proline stays trans on a straight run
    BOOST_TEST(peptide_side(3, 4) < 0);
    int cis_bonds = 0;
    for (int number = 1; number < 15; ++number)
        cis_bonds += peptide_side(number, number + 1) > 0;
    BOOST_TEST(cis_bonds == 0);
    // At the outside of the turn, CA-CB points to the right.
    BOOST_TEST(named_positions.at(7).at("CB").x >
               named_positions.at(7).at("CA").x);
    BOOST_TEST(named_positions.at(8).at("CB").x >
               named_positions.at(8).at("CA").x);
    const auto& proline = named_positions.at(3);
    const auto& preceding = named_positions.at(2);
    BOOST_TEST(std::abs(proline.at("CA").y - proline.at("N").y) < 1e-5);
    BOOST_TEST(proline.at("CA").x > proline.at("N").x);
    BOOST_TEST((proline.at("CB").y - proline.at("CA").y) *
                   (preceding.at("CB").y - preceding.at("CA").y) >
               0);
    for (const auto& [first, second] :
         {std::pair{"CA", "CB"}, {"CB", "CG"}, {"CG", "CD"}, {"CD", "N"}}) {
        const auto bond = proline.at(second) - proline.at(first);
        BOOST_TEST(std::abs(std::hypot(bond.x, bond.y) - bond_length) < 1e-5);
    }
}

BOOST_AUTO_TEST_CASE(TestMonomerLayoutKeepsUnrecognizedResidues)
{
    // An inline SMILES monomer has no database template but must retain
    // every atom with finite coordinates in the guided layout.
    auto monomer_mol = helm::helm_to_rdkit(
        "PEPTIDE1{A.[CC(O)[C@H](N[*:1])C(=O)[*:2]].G}$$$$V2.0");
    auto mol = toAtomistic(*monomer_mol);
    const auto num_atoms = mol->getNumAtoms();
    compute2DCoords(*mol);
    BOOST_TEST(mol->getNumAtoms() == num_atoms);
    BOOST_TEST(mol->getConformer().getNumAtoms() == num_atoms);
    for (unsigned int i = 0; i < num_atoms; ++i) {
        const auto p = mol->getConformer().getAtomPos(i);
        BOOST_TEST(std::isfinite(p.x));
        BOOST_TEST(std::isfinite(p.y));
    }
}

BOOST_AUTO_TEST_CASE(TestTerminalAmineFollowsPeptideBackbone)
{
    // Check an [am] cap after both a short straight chain and a turned
    // chain; its attachment must keep a normal peptide bond length.
    for (const auto& helm_string :
         {"PEPTIDE1{A.M.[am]}$$$$",
          "PEPTIDE1{[dR].P.[dK].[dP].Q.Q.F.[dF].G.L.M.[am]}$$$$"}) {
        auto monomer_mol = helm::helm_to_rdkit(helm_string);
        auto atomistic_mol = toAtomistic(*monomer_mol);
        compute2DCoords(*atomistic_mol);
        const auto& conf = atomistic_mol->getConformer();
        std::optional<RDGeom::Point3D> carbonyl, terminal_nitrogen;
        for (const auto* atom : atomistic_mol->atoms()) {
            const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                atom->getMonomerInfo());
            if (!info)
                continue;
            auto name = info->getName();
            name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
            if (info->getResidueNumber() ==
                    static_cast<int>(monomer_mol->getNumAtoms()) - 1 &&
                name == "C")
                carbonyl = conf.getAtomPos(atom->getIdx());
            if (info->getResidueNumber() ==
                    static_cast<int>(monomer_mol->getNumAtoms()) &&
                name == "N")
                terminal_nitrogen = conf.getAtomPos(atom->getIdx());
        }
        BOOST_REQUIRE(carbonyl.has_value());
        BOOST_REQUIRE(terminal_nitrogen.has_value());
        const auto bond = *terminal_nitrogen - *carbonyl;
        BOOST_TEST(std::abs(std::hypot(bond.x, bond.y) - RDDepict::BOND_LEN) <
                   1e-5);
    }
}

std::map<int, std::map<std::string, RDGeom::Point3D>>
residue_atom_positions(const RDKit::ROMol& mol)
{
    std::map<int, std::map<std::string, RDGeom::Point3D>> positions;
    for (const auto* atom : mol.atoms()) {
        const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        if (!info)
            continue;
        auto name = info->getName();
        name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
        positions[info->getResidueNumber()][name] =
            mol.getConformer().getAtomPos(atom->getIdx());
    }
    return positions;
}

BOOST_AUTO_TEST_CASE(TestAcetylCapAndTerminalCysteineKeepBackboneAngles)
{
    // The acetyl cap and terminal disulfide-linked Cys must keep their
    // peptide angles while side-chain selection shortens the S-S bridge.
    const auto helm_string =
        R"HELM(PEPTIDE1{[ac].[*N[C@@H](C(=O)*)S* |$_R1;;;;;_R2;;_R3$|].F.F.[dW].K.T.F.C}$PEPTIDE1,PEPTIDE1,2:R3-9:R3$$$)HELM";
    auto monomers = helm::helm_to_rdkit(helm_string);
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    const auto cosine = [](const RDGeom::Point3D& a, const RDGeom::Point3D& b) {
        return (a.x * b.x + a.y * b.y) /
               (std::hypot(a.x, a.y) * std::hypot(b.x, b.y));
    };
    const auto& cap = positions.at(1);
    const auto& first = positions.at(2);
    const auto cap_bond = cap.at("C") - first.at("N");
    const auto first_n_ca = first.at("CA") - first.at("N");
    BOOST_TEST(std::abs(cosine(cap_bond, first_n_ca) + 0.5) < 1e-5);
    BOOST_TEST(std::abs(std::hypot(cap_bond.x, cap_bond.y) -
                        RDDepict::BOND_LEN) < 1e-5);
    // With this starting phase, acetyl O points below the backbone and the
    // methyl occupies the opposite side of its carbonyl.
    BOOST_TEST(cap.at("O").y < cap.at("C").y);

    const auto& terminal = positions.at(9);
    const auto ca_n = terminal.at("N") - terminal.at("CA");
    const auto ca_c = terminal.at("C") - terminal.at("CA");
    // Without a following peptide bond, C9 can use a regular 120-degree
    // N-CA-C angle instead of opening the terminal residue into linear cis.
    BOOST_TEST(std::abs(cosine(ca_n, ca_c) + 0.5) < 1e-5);
    bool checked_disulfide = false;
    for (const auto* bond : mol->bonds()) {
        if (bond->getBeginAtom()->getAtomicNum() != 16 ||
            bond->getEndAtom()->getAtomicNum() != 16)
            continue;
        const auto disulfide =
            mol->getConformer().getAtomPos(bond->getBeginAtomIdx()) -
            mol->getConformer().getAtomPos(bond->getEndAtomIdx());
        // Joint side-chain selection shortens this bridge while retaining
        // the already checked 120-degree backbone geometry.
        BOOST_TEST(std::hypot(disulfide.x, disulfide.y) <
                   RDDepict::BOND_LEN * 4.75);
        checked_disulfide = true;
    }
    BOOST_TEST(checked_disulfide);
}

BOOST_AUTO_TEST_CASE(TestAcetylCapGeometryWithProline)
{
    // A cis N-CA proline changes the cap's C-N direction. Keep the O-C-N
    // angle trigonal and the methyl opposite N across the C=O bond.
    for (const auto* helm :
         {"PEPTIDE1{[ac].P.A}$$$$V2.0", "PEPTIDE1{[ac].A.P}$$$$V2.0"}) {
        auto monomers = helm::helm_to_rdkit(helm);
        auto mol = toAtomistic(*monomers);
        compute2DCoords(*mol);
        const auto positions = residue_atom_positions(*mol);
        const auto& cap = positions.at(1);
        const auto& first = positions.at(2);
        const auto c_o = cap.at("O") - cap.at("C");
        const auto c_n = first.at("N") - cap.at("C");
        const auto cosine =
            (c_o.x * c_n.x + c_o.y * c_n.y) /
            (std::hypot(c_o.x, c_o.y) * std::hypot(c_n.x, c_n.y));
        BOOST_TEST_CONTEXT(helm)
        {
            BOOST_TEST(std::abs(cosine + 0.5) < 1e-5);
        }
        bool checked_methyl = false;
        for (const auto* carbonyl : mol->atoms()) {
            const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                carbonyl->getMonomerInfo());
            if (!info || info->getResidueNumber() != 1) {
                continue;
            }
            auto name = info->getName();
            name.erase(std::remove(name.begin(), name.end(), ' '), name.end());
            if (name != "C") {
                continue;
            }
            for (const auto* methyl : mol->atomNeighbors(carbonyl)) {
                const auto* methyl_info =
                    dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                        methyl->getMonomerInfo());
                if (methyl->getAtomicNum() != 6 || !methyl_info ||
                    methyl_info->getResidueNumber() != 1) {
                    continue;
                }
                const auto c_methyl =
                    mol->getConformer().getAtomPos(methyl->getIdx()) -
                    cap.at("C");
                const auto cross = [&](const RDGeom::Point3D& v) {
                    return c_o.x * v.y - c_o.y * v.x;
                };
                BOOST_TEST(cross(c_methyl) * cross(c_n) < 0);
                checked_methyl = true;
            }
        }
        BOOST_TEST(checked_methyl);
    }
}

BOOST_AUTO_TEST_CASE(TestInlinePeptideBackboneUsesR1AndR2)
{
    // Inline SMILES residues with R1 and R2 belong to the peptide path;
    // their named backbone atoms and neighboring bond lengths must survive.
    const auto helm_string =
        R"HELM(PEPTIDE1{[ac].[*N[C@@H](C(=O)*)C(C)(C)S* |$_R1;;;;;_R2;;;;;_R3$|].[*N[C@@H](C(=O)*)[C@@H](C)CC |$_R1;;;;;_R2;;;;$|].P.[dR].G.[dD].[*N[C@H](C(=O)*)Cc1ccc(cc1)OC |$_R1;;;;;_R2;;;;;;;;;$|].[dR].[dC].[am]}$PEPTIDE1,PEPTIDE1,2:R3-10:R3$$$)HELM";
    auto monomers = helm::helm_to_rdkit(helm_string);
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    for (int number : {2, 3, 8}) {
        for (const auto* name : {"N", "CA", "C", "O", "CB"})
            BOOST_TEST(positions.at(number).contains(name));
    }
    BOOST_TEST(positions.at(2).at("C").x > positions.at(2).at("N").x);
    BOOST_TEST(positions.at(3).at("C").x > positions.at(3).at("N").x);
    BOOST_TEST(positions.at(8).at("C").x < positions.at(8).at("N").x);

    for (const auto& [first, second] :
         {std::pair{1, 2}, std::pair{2, 3}, std::pair{3, 4}, std::pair{7, 8},
          std::pair{8, 9}}) {
        const auto bond =
            positions.at(second).at("N") - positions.at(first).at("C");
        BOOST_TEST(std::abs(std::hypot(bond.x, bond.y) - RDDepict::BOND_LEN) <
                   1e-5);
    }
}

BOOST_AUTO_TEST_CASE(TestTurnStartingAtSecondResidue)
{
    // With only one link before the bend, the chain must still find the
    // returning arm and maintain all peptide bond lengths.
    const auto helm_string =
        R"HELM(PEPTIDE1{C.[*N[C@H](C(=O)*)C(c1ccccc1)c2ccccc2 |$_R1;;;;;_R2;;;;;;;;;;;;;$|].F.V.N.C.P.R.G.[am]}$PEPTIDE1,PEPTIDE1,1:R3-6:R3$$$)HELM";
    auto monomers = helm::helm_to_rdkit(helm_string);
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    for (int number = 1; number <= 9; ++number) {
        const auto peptide =
            positions.at(number + 1).at("N") - positions.at(number).at("C");
        BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) -
                            RDDepict::BOND_LEN) < 1e-5);
    }
    // The monomer layout turns immediately after C1. The lower arm must
    // return leftward, even though only one link precedes the bend.
    BOOST_TEST(positions.at(3).at("CA").x > positions.at(1).at("CA").x);
    BOOST_TEST(positions.at(6).at("CA").x < positions.at(4).at("CA").x);
    BOOST_TEST(positions.at(8).at("CA").x < positions.at(6).at("CA").x);
}

BOOST_AUTO_TEST_CASE(TestUnnamedInlinePeptideCarbonylFollowsBackbone)
{
    // Two chains begin with an unnamed SMILES carbonyl. Both attachments
    // must have a normal C-N length and keep the carbon branch behind it.
    const auto helm_string =
        R"HELM(PEPTIDE1{[[C@H](C(=O)*)CCCN |$;;;_R2;;;;$|].C.[dV]}|PEPTIDE2{[[C@H](C(=O)*)CCCN |$;;;_R2;;;;$|].C.[dV]}$PEPTIDE1,PEPTIDE2,2:R3-2:R3$$$)HELM";
    auto monomers = helm::helm_to_rdkit(helm_string);
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto& conf = mol->getConformer();
    int attached_custom_residues = 0;
    for (const auto* bond : mol->bonds()) {
        const auto* first = bond->getBeginAtom();
        const auto* second = bond->getEndAtom();
        const auto* first_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            first->getMonomerInfo());
        const auto* second_info =
            dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                second->getMonomerInfo());
        if (!first_info || !second_info ||
            first_info->getChainId() != second_info->getChainId())
            continue;
        const RDKit::Atom* carbonyl = nullptr;
        const RDKit::Atom* nitrogen = nullptr;
        if (first_info->getResidueNumber() == 1 &&
            second_info->getResidueNumber() == 2) {
            carbonyl = first;
            nitrogen = second;
        } else if (second_info->getResidueNumber() == 1 &&
                   first_info->getResidueNumber() == 2) {
            carbonyl = second;
            nitrogen = first;
        } else {
            continue;
        }
        BOOST_REQUIRE(carbonyl->getAtomicNum() == 6);
        BOOST_REQUIRE(nitrogen->getAtomicNum() == 7);
        const auto peptide = conf.getAtomPos(nitrogen->getIdx()) -
                             conf.getAtomPos(carbonyl->getIdx());
        BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) -
                            RDDepict::BOND_LEN) < 1e-5);
        const RDKit::Atom* preceding_carbon = nullptr;
        for (const auto* neighbor : mol->atomNeighbors(carbonyl)) {
            const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                neighbor->getMonomerInfo());
            if (neighbor->getAtomicNum() == 6 && info &&
                info->getChainId() == first_info->getChainId() &&
                info->getResidueNumber() == 1)
                preceding_carbon = neighbor;
        }
        BOOST_REQUIRE(preceding_carbon);
        const auto branch = conf.getAtomPos(preceding_carbon->getIdx()) -
                            conf.getAtomPos(carbonyl->getIdx());
        BOOST_TEST(branch.x * peptide.x + branch.y * peptide.y < 0);
        ++attached_custom_residues;
    }
    BOOST_TEST(attached_custom_residues == 2);
}

BOOST_AUTO_TEST_CASE(TestTerminalInlinePeptideWithoutR2)
{
    // A final SMILES residue lacking R2 remains part of the backbone,
    // including its carbonyl and the incoming peptide bond.
    const auto helm_string =
        R"HELM(PEPTIDE1{[ac].C.F.F.[dW].K.T.F.[*N[C@H](C(=O)O)S* |$_R1;;;;;;;_R3$|]}$PEPTIDE1,PEPTIDE1,2:R3-9:R3$$$)HELM";
    auto monomers = helm::helm_to_rdkit(helm_string);
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    for (const auto* name : {"N", "CA", "C", "O", "CB"})
        BOOST_TEST(positions.at(9).contains(name));
    const auto peptide = positions.at(9).at("N") - positions.at(8).at("C");
    BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) - RDDepict::BOND_LEN) <
               1e-5);
    // The unnamed acid oxygen must move with the terminal carbonyl carbon.
    bool checked_acid_oxygen = false;
    for (const auto* bond : mol->bonds()) {
        for (const auto& [carbonyl, oxygen] :
             {std::pair{bond->getBeginAtom(), bond->getEndAtom()},
              std::pair{bond->getEndAtom(), bond->getBeginAtom()}}) {
            const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                oxygen->getMonomerInfo());
            if (carbonyl->getAtomicNum() != 6 || oxygen->getAtomicNum() != 8 ||
                bond->getBondType() != RDKit::Bond::SINGLE || !info ||
                info->getResidueNumber() != 9 || !info->getName().empty()) {
                continue;
            }
            const auto delta =
                mol->getConformer().getAtomPos(oxygen->getIdx()) -
                mol->getConformer().getAtomPos(carbonyl->getIdx());
            BOOST_TEST(std::abs(std::hypot(delta.x, delta.y) -
                                RDDepict::BOND_LEN * 0.7) < 1e-5);
            checked_acid_oxygen = true;
        }
    }
    BOOST_TEST(checked_acid_oxygen);
}

BOOST_AUTO_TEST_CASE(TestTerminalInlinePeptideAmideFollowsTurn)
{
    // A final SMILES residue can contain its own amide N and lack R2;
    // check its backbone, return arm, and terminal C-N length.
    const auto helm_string =
        R"HELM(PEPTIDE1{[ac].[dC].G.[*N(C)[C@@H](C(=O)*)CCCN=C(N)N |$_R1;;;;;;_R2;;;;;;;$|].G.[dD].[*N[C@@H](C(=O)N)C(C)(C)S* |$_R1;;;;;;;;;;_R3$|]}$PEPTIDE1,PEPTIDE1,2:R3-7:R3$$$)HELM";
    auto monomers = helm::helm_to_rdkit(helm_string);
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);

    // The terminal amide has no R2. Its N-CA-C(=O) backbone still belongs
    // to the peptide chain, so every peptide bond keeps a normal length.
    for (int number : {4, 7})
        for (const auto* name : {"N", "CA", "C", "O"})
            BOOST_TEST(positions.at(number).contains(name));
    for (int number = 1; number < 7; ++number) {
        const auto peptide =
            positions.at(number + 1).at("N") - positions.at(number).at("C");
        BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) -
                            RDDepict::BOND_LEN) < 1e-5);
    }
    BOOST_TEST(positions.at(7).at("CA").x < positions.at(4).at("CA").x);
    bool checked_amide = false;
    const auto& conf = mol->getConformer();
    for (const auto* bond : mol->bonds()) {
        const auto* first = bond->getBeginAtom();
        const auto* second = bond->getEndAtom();
        for (const auto& [carbonyl, amide] :
             {std::pair{first, second}, std::pair{second, first}}) {
            const auto* info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                amide->getMonomerInfo());
            if (carbonyl->getAtomicNum() != 6 || amide->getAtomicNum() != 7 ||
                !info || info->getResidueNumber() != 7 ||
                !info->getName().empty())
                continue;
            const auto delta = conf.getAtomPos(amide->getIdx()) -
                               conf.getAtomPos(carbonyl->getIdx());
            BOOST_TEST(std::abs(std::hypot(delta.x, delta.y) -
                                RDDepict::BOND_LEN * 0.7) < 1e-5);
            checked_amide = true;
        }
    }
    BOOST_TEST(checked_amide);
}

BOOST_AUTO_TEST_CASE(TestDisulfideCysteineFacesItsPartner)
{
    // Starting phase should face both linked Cys side chains toward the
    // bridge without opening an extra residue on the straight arm.
    // Flipping the starting phase faces both sides of this bridge without
    // opening C7 into LINEAR_CIS_N_CA.
    auto monomers =
        helm::helm_to_rdkit("PEPTIDE1{R.G.G.R.L.C.C.C.R.R.R.F.C.V.C.V.G.R}"
                            "$PEPTIDE1,PEPTIDE1,7:R3-13:R3$$$");
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    const auto faces_partner = [&](int number, int partner) {
        const auto& residue = positions.at(number);
        const auto backbone = residue.at("C") - residue.at("N");
        const auto sidechain = residue.at("CB") - residue.at("CA");
        const auto toward_partner =
            positions.at(partner).at("CA") - residue.at("CA");
        const auto cross = [&](const RDGeom::Point3D& direction) {
            return backbone.x * direction.y - backbone.y * direction.x;
        };
        return cross(sidechain) * cross(toward_partner) > 0;
    };
    BOOST_TEST(faces_partner(7, 13));
    BOOST_TEST(faces_partner(13, 7));
    BOOST_TEST(std::abs(positions.at(7).at("CA").y -
                        positions.at(7).at("N").y) > RDDepict::BOND_LEN / 4);
    BOOST_TEST(std::abs(positions.at(8).at("CA").y -
                        positions.at(8).at("N").y) > RDDepict::BOND_LEN / 4);
    const auto side = [&](int number) {
        return positions.at(number).at("CB").y -
               positions.at(number).at("CA").y;
    };
    BOOST_TEST(side(6) * side(7) < 0); // C7 stays ALL_TRANS
    BOOST_TEST(side(7) * side(8) < 0); // C8 remains ALL_TRANS
    const auto upper = positions.at(4).at("C") - positions.at(4).at("N");
    const auto lower = positions.at(16).at("C") - positions.at(16).at("N");
    BOOST_TEST(std::abs(upper.x * lower.y - upper.y * lower.x) < 1e-5);
    BOOST_TEST(upper.x * lower.x + upper.y * lower.y < 0);
    for (int number = 6; number <= 8; ++number) {
        const auto peptide =
            positions.at(number + 1).at("N") - positions.at(number).at("C");
        BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) -
                            RDDepict::BOND_LEN) < 1e-5);
    }
}

BOOST_AUTO_TEST_CASE(TestDProlineInMixedPeptide)
{
    // dProline in a mixed sequence needs the same straight-chain ring
    // geometry and preceding-side-chain direction as ordinary proline.
    auto monomers = helm::helm_to_rdkit(
        "PEPTIDE1{[dR].P.[dK].[dP].Q.Q.F.[dF].G.L.M.[am]}$$$$");
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    const auto& preceding = positions.at(3);
    const auto& dproline = positions.at(4);
    BOOST_TEST(std::abs(dproline.at("N").y - dproline.at("CA").y) < 1e-5);
    BOOST_TEST(dproline.at("CA").x > dproline.at("N").x);
    BOOST_TEST((preceding.at("CB").y - preceding.at("CA").y) *
                   (dproline.at("CB").y - dproline.at("CA").y) >
               0);
}

BOOST_AUTO_TEST_CASE(TestFourResidueTurnSpreadsAngles)
{
    // A four-residue D-peptide bend must preserve each turn direction while
    // spreading its angle correction and keeping peptide links trans.
    auto monomers =
        helm::helm_to_rdkit("PEPTIDE1{[dK].[dF].[dA].[dK].[dF].[dA].[dK].[dK].["
                            "dF].[dA].[dK].[dF].[dA].[dK]}$$$$");
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    BOOST_REQUIRE(positions.size() == 14);
    const auto tangent = [&](int number) {
        return positions.at(number).at("C") - positions.at(number).at("N");
    };
    const auto signed_turn = [&](int number) {
        const auto a = tangent(number), b = tangent(number + 1);
        return std::atan2(a.x * b.y - a.y * b.x, a.x * b.x + a.y * b.y);
    };
    // The four residues around the bend retain the original direction of
    // each backbone turn. Opening the bend changes each local angle only
    // slightly, so the two straight arms become antiparallel.
    for (int number = 5; number <= 8; ++number) {
        BOOST_TEST(signed_turn(number) < 0);
        BOOST_TEST(std::abs(signed_turn(number)) < std::numbers::pi / 3);
    }
    BOOST_TEST((positions.at(6).at("C") - positions.at(6).at("CA")).x > 0);
    const auto upper = tangent(3), lower = tangent(12);
    BOOST_TEST(std::abs(upper.x * lower.y - upper.y * lower.x) < 1e-5);
    BOOST_TEST(upper.x * lower.x + upper.y * lower.y < 0);

    std::vector<int> cis_links;
    for (int number = 1; number < 14; ++number) {
        const auto& current = positions.at(number);
        const auto& next = positions.at(number + 1);
        const auto peptide = next.at("N") - current.at("C");
        BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) -
                            RDDepict::BOND_LEN) < 1e-5);
        const auto before_alpha = current.at("CA") - current.at("C");
        const auto after_alpha = next.at("CA") - next.at("N");
        const auto cross = [&](const RDGeom::Point3D& direction) {
            return peptide.x * direction.y - peptide.y * direction.x;
        };
        if (cross(before_alpha) * cross(after_alpha) > 0)
            cis_links.push_back(number);
    }
    BOOST_TEST(cis_links.empty());
}

BOOST_AUTO_TEST_CASE(TestTurnsKeepArmsAntiparallel)
{
    // Both an early proline and S-S bridges can complicate turn detection.
    // In each chain, the straight arms must remain antiparallel and every
    // peptide bond must retain its normal length.
    struct TurnCase {
        const char* helm;
        int upper_residue;
        int lower_residue;
        int last_peptide_residue;
    };
    const TurnCase cases[] = {
        {"PEPTIDE1{L.P.N.Y.N.W.N.S.F.G.I.R.F.[am]}$$$$", 5, 10, 13},
        {"PEPTIDE1{R.G.G.R.L.C.Y.C.R.R.R.F.C.V.C.V.G.R}"
         "$PEPTIDE1,PEPTIDE1,8:R3-13:R3|PEPTIDE1,PEPTIDE1,6:R3-15:R3$$$",
         4, 16, 17},
    };
    for (const auto& test_case : cases) {
        BOOST_TEST_CONTEXT(test_case.helm)
        {
            auto monomers = helm::helm_to_rdkit(test_case.helm);
            auto mol = toAtomistic(*monomers);
            compute2DCoords(*mol);
            const auto positions = residue_atom_positions(*mol);
            const auto direction = [&](int number) {
                const auto& residue = positions.at(number);
                return residue.at("C") - residue.at("N");
            };
            const auto upper = direction(test_case.upper_residue);
            const auto lower = direction(test_case.lower_residue);
            BOOST_TEST(std::abs(upper.x * lower.y - upper.y * lower.x) < 1e-5);
            BOOST_TEST(upper.x * lower.x + upper.y * lower.y < 0);
            for (int number = 1; number <= test_case.last_peptide_residue;
                 ++number) {
                const auto peptide = positions.at(number + 1).at("N") -
                                     positions.at(number).at("C");
                BOOST_TEST(std::abs(std::hypot(peptide.x, peptide.y) -
                                    RDDepict::BOND_LEN) < 1e-5);
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(TestProlineInTurnPrefersFlippedPhase)
{
    // The opposite starting phase avoids an unnecessary linear cis residue
    // before a proline turn while preserving the returning arm.
    auto monomers =
        helm::helm_to_rdkit("PEPTIDE1{[ac].S.N.W.K.W.W.P.G.I.F.D.[am]}$$$$");
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    const auto direction = [&](int number) {
        const auto& residue = positions.at(number);
        return residue.at("C") - residue.at("N");
    };
    const auto& lysine = positions.at(5);
    // The opposite starting phase keeps K trans, avoiding the linear cis
    // correction before the W-W-P bend. The bend still turns at N-CA.
    BOOST_TEST(std::abs(lysine.at("CA").y - lysine.at("N").y) >
               RDDepict::BOND_LEN / 4);
    BOOST_TEST(std::abs(lysine.at("C").y - lysine.at("CA").y) >
               RDDepict::BOND_LEN / 4);
    const auto upper = direction(4), lower = direction(10);
    BOOST_TEST(std::abs(upper.x * lower.y - upper.y * lower.x) < 1e-5);
    BOOST_TEST(upper.x * lower.x + upper.y * lower.y < 0);
    for (int number = 6; number <= 7; ++number) {
        const auto a = direction(number), b = direction(number + 1);
        const auto turn =
            std::atan2(a.x * b.y - a.y * b.x, a.x * b.x + a.y * b.y);
        BOOST_TEST(std::abs(turn + std::numbers::pi / 3) < 1e-5);
    }
}

BOOST_AUTO_TEST_CASE(TestThreeResidueTurnNeedsNoAngleCorrection)
{
    // An eleven-alanine chain makes a three-residue bend that should keep
    // its regular 60-degree turns without angle spreading.
    auto monomers = helm::helm_to_rdkit("PEPTIDE1{A.A.A.A.A.A.A.A.A.A.A}$$$$");
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    const auto direction = [&](int number) {
        const auto& residue = positions.at(number);
        return residue.at("C") - residue.at("N");
    };
    const auto upper = direction(4), lower = direction(8);
    BOOST_TEST(std::abs(upper.x * lower.y - upper.y * lower.x) < 1e-5);
    BOOST_TEST(upper.x * lower.x + upper.y * lower.y < 0);
    for (int number = 5; number <= 7; ++number) {
        const auto a = direction(number), b = direction(number + 1);
        const auto turn =
            std::atan2(a.x * b.y - a.y * b.x, a.x * b.x + a.y * b.y);
        BOOST_TEST(std::abs(turn + std::numbers::pi / 3) < 1e-5);
    }
}

BOOST_AUTO_TEST_CASE(TestDisulfidePrefersFlippedStartingPhase)
{
    // The flipped start lets the first dC and a later Cys face their S-S
    // partners without changing the straight backbone geometry.
    auto monomers = helm::helm_to_rdkit(
        "PEPTIDE1{[dC].[dS].C.[dS].S.[dL].M.[dD].K.[dE].C.[dV].Y.[dF].[dC].H.A."
        "D.I.I.W}"
        "$PEPTIDE1,PEPTIDE1,3:R3-11:R3|PEPTIDE1,PEPTIDE1,1:R3-15:R3$$$");
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    const auto& first = positions.at(1);
    const auto& cysteine = positions.at(3);
    const auto& partner = positions.at(11);
    BOOST_TEST(std::abs(first.at("CA").y - first.at("N").y) >
               RDDepict::BOND_LEN / 4);
    BOOST_TEST(first.at("CB").y < first.at("CA").y);
    BOOST_TEST(std::abs(cysteine.at("CA").y - cysteine.at("N").y) >
               RDDepict::BOND_LEN / 4);
    BOOST_TEST(cysteine.at("SG").x > cysteine.at("CA").x);
    const auto backbone = cysteine.at("C") - cysteine.at("N");
    const auto sidechain = cysteine.at("SG") - cysteine.at("CA");
    const auto toward_partner = partner.at("CA") - cysteine.at("CA");
    const auto cross = [&](const RDGeom::Point3D& direction) {
        return backbone.x * direction.y - backbone.y * direction.x;
    };
    BOOST_TEST(cross(sidechain) * cross(toward_partner) > 0);
}

BOOST_AUTO_TEST_CASE(TestTurnKeepsArmsApart)
{
    // In this Arg-rich turn, angle spreading should separate the two arms
    // while retaining the direction of the local backbone turns.
    auto monomers =
        helm::helm_to_rdkit("PEPTIDE1{T.P.Q.R.A.R.R.R.K.K.R.F}$$$$");
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    const auto direction = [&](int number) {
        const auto& residue = positions.at(number);
        return residue.at("C") - residue.at("N");
    };
    const auto upper = direction(4), lower = direction(9);
    BOOST_TEST(std::abs(upper.x * lower.y - upper.y * lower.x) < 1e-5);
    BOOST_TEST(upper.x * lower.x + upper.y * lower.y < 0);
    // A5, R6, R7 and R8 form the bend. Both arms must remain separated,
    // with the backbone continuing to turn in the same direction.
    BOOST_TEST(positions.at(9).at("N").y <
               positions.at(4).at("N").y - 2 * RDDepict::BOND_LEN);
    for (int number = 6; number <= 7; ++number) {
        const auto a = direction(number), b = direction(number + 1);
        const auto turn =
            std::atan2(a.x * b.y - a.y * b.x, a.x * b.x + a.y * b.y);
        BOOST_TEST(turn < 0);
        BOOST_TEST(std::abs(turn + std::numbers::pi / 4) <
                   std::numbers::pi / 12);
    }
    const auto last_loop = direction(8), first_lower = direction(9);
    BOOST_TEST(std::abs(last_loop.x * first_lower.y -
                        last_loop.y * first_lower.x) < 1e-5);
}

BOOST_AUTO_TEST_CASE(TestTwoTurnsKeepMiddleArmStraight)
{
    // Two bends in one chain should leave the intervening arm straight
    // rather than creating an extra turn in its middle.
    auto monomers = helm::helm_to_rdkit("PEPTIDE1{H.S.D.G.I.A.T.D.S.Y.S.R.Y.R."
                                        "K.Q.M.A.V.K.K.Y.L.A.A.V.L.[am]}$$$$");
    auto mol = toAtomistic(*monomers);
    compute2DCoords(*mol);
    const auto positions = residue_atom_positions(*mol);
    const auto direction = [&](int number) {
        const auto& residue = positions.at(number);
        return residue.at("C") - residue.at("N");
    };
    const auto upper = direction(5);
    const auto middle = direction(14);
    const auto lower = direction(24);
    BOOST_TEST(std::abs(upper.x * middle.y - upper.y * middle.x) < 1e-5);
    BOOST_TEST(std::abs(middle.x * lower.y - middle.y * lower.x) < 1e-5);
    BOOST_TEST(upper.x * middle.x + upper.y * middle.y < 0);
    BOOST_TEST(middle.x * lower.x + middle.y * lower.y < 0);
    for (int number = 12; number <= 17; ++number) {
        const auto current = direction(number);
        BOOST_TEST(std::abs(current.x * middle.y - current.y * middle.x) <
                   1e-5);
        BOOST_TEST(current.x * middle.x + current.y * middle.y > 0);
    }
    BOOST_TEST(positions.at(5).at("N").y >
               positions.at(14).at("N").y + 2 * RDDepict::BOND_LEN);
    BOOST_TEST(positions.at(14).at("N").y >
               positions.at(24).at("N").y + 2 * RDDepict::BOND_LEN);
}

bool has_atom_map_numbers(const RDKit::ROMol& mol)
{
    for (const auto* atom : mol.atoms()) {
        if (atom->hasProp(RDKit::common_properties::molAtomMapNumber)) {
            return true;
        }
    }
    return false;
};

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

BOOST_AUTO_TEST_CASE(TestCyclicPolymerWithMissingMonomerInAnotherPolymer)
{
    auto monomer_mol = to_rdkit(
        "PEPTIDE1{A.A.A.A}|PEPTIDE2{A.A}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0");
    auto atomistic_mol = toAtomistic(*monomer_mol);
    atomistic_mol->setProp("SEQRES",
                           R"({"A":["A","A","A","A"],"B":["A","A","A"]})");

    auto roundtrip_monomer_mol = toMonomeric(*atomistic_mol);

    BOOST_TEST(
        to_string(*roundtrip_monomer_mol, Format::HELM) ==
        "PEPTIDE1{A.A.A.A}|PEPTIDE2{A.A.A}$PEPTIDE1,PEPTIDE1,4:R2-1:R1$$$V2.0");
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
         "S.G.V.P.V.P.G.V.W.L.V}|PEPTIDE2{[CCCCCC[C@H]1C[C@H]1CCCCCCCCCC(O)O.[*"
         ":1]]}|PEPTIDE3{V.S.G.T.W.Y.L.K.A.M.T.V.D.V.G.A.L.R.C.L.A.G.S.V.I.P.T."
         "T.L.T.T}|PEPTIDE4{N.L.E.A.K.V.T}|PEPTIDE5{E.V.K.A.V.L.S.K.T.D.E.P.G."
         "I.Y.T.A.I.G.G.I.H.V.A.K.I.G.R.S}|PEPTIDE6{D.H.Y.I.F.Y.S.E.G.C.L.S.G."
         "V.P.V.P.G.V.W.L.V}|PEPTIDE7{[CCCCCC[C@H]1C[C@H]1CCCCCCCCCC(O)O.[*:1]]"
         "}|PEPTIDE8{E.V.V.K.F.M.D.V.Y.Q.R.S.Y.C.H.P.I.E.T.L.V.D.I.F.Q.E.Y.P.D."
         "E.I.E.Y.I.F.K.P.S.C.V.P.L.M.R.C.G.G.C.C.N.D.E.G.L.E.C.V.P.T.E.E.S.N."
         "I.T.M.Q.I.M.R.I.K.P.H.Q.G.Q.H.I.G.E.M.S.F.L.Q.H.N.K.C.E.C.R.P.K.K}|"
         "PEPTIDE9{[OS(O)(O)O.[*:1]]}|PEPTIDE10{[OS(O)(O)O.[*:1]]}|PEPTIDE11{H."
         "E.V.V.K.F.M.D.V.Y.Q.R.S.Y.C.H.P.I.E.T.L.V.D.I.F.Q.E.Y.P.D.E.I.E.Y.I."
         "F.K.P.S.C.V.P.L.M.R.C.G.G.C.C.N.D.E.G.L.E.C.V.P.T.E.E.S.N.I.T.M.Q.I."
         "M.R.I.K.P.H.Q.G.Q.H.I.G.E.M.S.F.L.Q.H.N.K.C.E.C.R.P}|PEPTIDE12{[OS(O)"
         "(O)O.[*:1]]}|PEPTIDE13{[CC(O)O.[*:1]]}$PEPTIDE1,PEPTIDE1,24:R3-97:R3|"
         "PEPTIDE8,PEPTIDE11,39:R3-49:R3|PEPTIDE8,PEPTIDE8,45:R3-90:R3|"
         "PEPTIDE8,PEPTIDE11,48:R3-40:R3|PEPTIDE8,PEPTIDE8,49:R3-92:R3|"
         "PEPTIDE11,PEPTIDE11,46:R3-91:R3|PEPTIDE11,PEPTIDE11,50:R3-93:R3|"
         "PEPTIDE3,PEPTIDE6,19:R3-10:R3$$$V2.0"},
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

BOOST_AUTO_TEST_CASE(TestHelmImportWithCustomMonomer)
{
    // Register custom monomers matching the pattern used in the wasm tests:
    // multi-character symbols, realistic SMILES with attachment points.
    constexpr std::string_view custom_json = R"JSON([
        {
            "symbol": "Sar",
            "polymer_type": "PEPTIDE",
            "natural_analog": "G",
            "smiles": "[H:3]N(C)C([H:1])C(=O)[OH:2]",
            "core_smiles": "CN(CC(O)=O)[H]",
            "name": "Sarcosine",
            "monomer_type": "Backbone",
            "author": "test"
        },
        {
            "symbol": "TestMon",
            "polymer_type": "PEPTIDE",
            "natural_analog": "A",
            "smiles": "CC(C)(N[H:1])C(=O)[OH:2]",
            "core_smiles": "CC(C)(N)C=O",
            "name": "Test Monomer",
            "monomer_type": "Backbone",
            "author": "test"
        }
    ])JSON";

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(custom_json);

    // Test importing HELM strings that reference custom monomers.
    // These mirror the wasm tests that use sketcher_import_text with
    // autodetect format.
    const std::vector<std::string> helm_strings = {
        "PEPTIDE1{A.[Sar].G}$$$$V2.0",
        "PEPTIDE1{[TestMon].A}$$$$V2.0",
        "PEPTIDE1{A.[Sar].[TestMon].G}$$$$V2.0",
    };

    for (const auto& helm_str : helm_strings) {
        BOOST_TEST_CONTEXT("HELM: " << helm_str)
        {
            // Verify custom monomers pass validation
            auto validated_mol = helm::helm_to_rdkit(
                helm_str, /*do_throw=*/true, /*validate=*/true);
            BOOST_REQUIRE(validated_mol != nullptr);

            // Autodetect import (the path used by sketcher_import_text)
            auto monomer_mol = to_rdkit(helm_str, Format::AUTO_DETECT);
            BOOST_REQUIRE(monomer_mol != nullptr);

            // Convert to atomistic — requires DB lookup of custom monomer
            auto atomistic_mol = toAtomistic(*monomer_mol);
            BOOST_REQUIRE(atomistic_mol != nullptr);
            BOOST_TEST(atomistic_mol->getNumAtoms() > 0);

            // Roundtrip back to HELM
            auto roundtrip_monomer = toMonomeric(*atomistic_mol);
            BOOST_REQUIRE(roundtrip_monomer != nullptr);

            auto helm_result = to_string(*roundtrip_monomer, Format::HELM);
            BOOST_TEST(!helm_result.empty());
        }
    }
}

BOOST_AUTO_TEST_CASE(TestAtomisticConversionWithMultiAtomLeavingGroupSmiles)
{
    ResetMonomerDatabase reset_db;

    constexpr std::string_view custom_json = R"JSON([
        {
            "symbol": "AzA",
            "polymer_type": "PEPTIDE",
            "natural_analog": "X",
            "smiles": "[N-]=[N+]=[N:4]CC(=O)[OH:1]",
            "name": "Azidoacetic acid",
            "monomer_type": "Backbone",
            "author": "Example",
            "pdbcode": "AZA"
        }
    ])JSON";

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(custom_json);

    auto uncapped_monomer = to_rdkit("PEPTIDE1{A.[AzA]}$$$$V2.0", Format::HELM);
    BOOST_CHECK_THROW(toAtomistic(*uncapped_monomer), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestAtomisticConversionWithMultiAtomLeavingGroupAzide)
{
    ResetMonomerDatabase reset_db;

    constexpr std::string_view custom_json = R"JSON([
        {
            "symbol": "AzA",
            "polymer_type": "PEPTIDE",
            "natural_analog": "X",
            "smiles": "[N-]=[N+]=[N:4]CC(=O)[OH:1]",
            "core_smiles": "*CC(*)=O |$_R4;;;_R1;$|",
            "rgroups": "{\"R1\":\"*O |$_R1;$|\",\"R4\":\"*N=[N+]=[N-] |$_R4;;;$|\"}",
            "name": "Azidoacetic acid",
            "monomer_type": "Backbone",
            "author": "Example",
            "pdbcode": "AZA"
        }
    ])JSON";

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(custom_json);

    // No connection, azide group present but no atom map numbers
    auto uncapped_monomer = to_rdkit("PEPTIDE1{A.[AzA]}$$$$V2.0", Format::HELM);
    auto capped_atomistic = toAtomistic(*uncapped_monomer);
    std::unique_ptr<RDKit::RWMol> aza_query(RDKit::SmartsToMol("N=[N+]=[N-]"));
    auto capped_matches = RDKit::SubstructMatch(*capped_atomistic, *aza_query);
    BOOST_TEST(capped_matches.size() == 1);
    BOOST_TEST(!has_atom_map_numbers(*capped_atomistic));

    // Connection between AzA R4 and Alanine R1 carbon
    // No azide group and no atom map numbers
    auto linked_monomer = to_rdkit(
        "PEPTIDE1{A.[AzA]}$PEPTIDE1,PEPTIDE1,1:R1-2:R4$$$V2.0", Format::HELM);
    auto linked_atomistic = toAtomistic(*linked_monomer);
    auto linked_matches = RDKit::SubstructMatch(*linked_atomistic, *aza_query);
    BOOST_TEST(linked_matches.size() == 0);
    BOOST_TEST(!has_atom_map_numbers(*linked_atomistic));
}

BOOST_AUTO_TEST_CASE(TestAtomisticConversionWithMultiAtomLeavingGroupVinyl)
{
    ResetMonomerDatabase reset_db;

    constexpr std::string_view custom_json = R"JSON([
        {
            "symbol": "VnY",
            "polymer_type": "PEPTIDE",
            "natural_analog": "X",
            "smiles": "[C]=[C:4]CC(=O)[OH:1]",
            "core_smiles": "*CC(*)=O |$_R4;;;_R1;$|",
            "rgroups": "{\"R1\":\"*O |$_R1;$|\",\"R4\":\"*C=[C] |$_R4;;;$|\"}",
            "name": "Vinylacetic acid",
            "monomer_type": "Backbone",
            "author": "Example",
            "pdbcode": "VNY"
        }
    ])JSON";

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(custom_json);

    // No connection, vinyl group should be present but no atom map nums
    auto uncapped_monomer = to_rdkit("PEPTIDE1{A.[VnY]}$$$$V2.0", Format::HELM);
    auto capped_atomistic = toAtomistic(*uncapped_monomer);
    std::unique_ptr<RDKit::RWMol> vinyl_query(RDKit::SmilesToMol("C=C"));
    auto capped_matches =
        RDKit::SubstructMatch(*capped_atomistic, *vinyl_query);
    BOOST_TEST(capped_matches.size() == 1);
    BOOST_TEST(!has_atom_map_numbers(*capped_atomistic));

    // Connection between VnY R4 and Alanine R1 carbon
    // No vinyl group should be present and no atom map nums
    auto linked_monomer = to_rdkit(
        "PEPTIDE1{A.[VnY]}$PEPTIDE1,PEPTIDE1,1:R1-2:R4$$$V2.0", Format::HELM);
    auto linked_atomistic = toAtomistic(*linked_monomer);
    auto linked_matches =
        RDKit::SubstructMatch(*linked_atomistic, *vinyl_query);
    BOOST_TEST(linked_matches.size() == 0);
    BOOST_TEST(!has_atom_map_numbers(*linked_atomistic));
}

BOOST_AUTO_TEST_CASE(TestAtomisticConversionWithMultiAtomLeavingGroupAllyl)
{
    ResetMonomerDatabase reset_db;

    constexpr std::string_view custom_json = R"JSON([
        {
            "symbol": "AlY",
            "polymer_type": "PEPTIDE",
            "natural_analog": "X",
            "smiles": "[C]=[C]-[C:4]CC(=O)[OH:1]",
            "core_smiles": "*CC(*)=O |$_R4;;;_R1;$|",
            "rgroups": "{\"R1\":\"*O |$_R1;$|\",\"R4\":\"*C-[C]=[C] |$_R4;;;$|\"}",
            "name": "Allyl acid",
            "monomer_type": "Backbone",
            "author": "Example",
            "pdbcode": "ALY"
        }
    ])JSON";

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(custom_json);

    // No connection, allyl group should be present but no atom map nums
    auto uncapped_monomer = to_rdkit("PEPTIDE1{A.[AlY]}$$$$V2.0", Format::HELM);
    auto capped_atomistic = toAtomistic(*uncapped_monomer);
    std::unique_ptr<RDKit::RWMol> allyl_query(RDKit::SmilesToMol("C-C=C"));
    auto capped_matches =
        RDKit::SubstructMatch(*capped_atomistic, *allyl_query);
    BOOST_TEST(capped_matches.size() == 1);
    BOOST_TEST(!has_atom_map_numbers(*capped_atomistic));

    // Connection between AlY R4 and Alanine R1 carbon
    // No allyl group should be present and no atom map nums
    auto linked_monomer = to_rdkit(
        "PEPTIDE1{A.[AlY]}$PEPTIDE1,PEPTIDE1,1:R1-2:R4$$$V2.0", Format::HELM);
    auto linked_atomistic = toAtomistic(*linked_monomer);
    auto linked_matches =
        RDKit::SubstructMatch(*linked_atomistic, *allyl_query);
    BOOST_TEST(linked_matches.size() == 0);
    BOOST_TEST(!has_atom_map_numbers(*linked_atomistic));
}

BOOST_AUTO_TEST_CASE(TestAtomisticConversionWithMultiAtomLeavingGroupMultiMatch)
{
    ResetMonomerDatabase reset_db;

    constexpr std::string_view custom_json = R"JSON([
        {
            "symbol": "TriC",
            "polymer_type": "PEPTIDE",
            "natural_analog": "X",
            "smiles": "CC[CH2:4]C(=O)[OH:1]",
            "core_smiles": "*C(*)=O |$_R4;;;_R1;$|",
            "rgroups": "{\"R1\":\"*O |$_R1;$|\",\"R4\":\"*CCC |$_R4;;;$|\"}",
            "name": "Propane",
            "monomer_type": "Backbone",
            "author": "Example",
            "pdbcode": "TRC"
        }
    ])JSON";

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(custom_json);

    // No connection, TriC group present but no atom map numbers
    auto uncapped_monomer =
        to_rdkit("PEPTIDE1{A.[TriC]}$$$$V2.0", Format::HELM);
    auto capped_atomistic = toAtomistic(*uncapped_monomer);
    std::unique_ptr<RDKit::RWMol> tri_query(RDKit::SmartsToMol("CCC"));
    auto capped_matches = RDKit::SubstructMatch(*capped_atomistic, *tri_query);
    // 5 matches between alanine and tric (7 carbons in a chain)
    BOOST_TEST(capped_matches.size() == 5);
    BOOST_TEST(!has_atom_map_numbers(*capped_atomistic));

    // Connection between Tri R4 and Alanine R1 carbon
    // No TriC group and no atom map numbers
    auto linked_monomer = to_rdkit(
        "PEPTIDE1{A.[TriC]}$PEPTIDE1,PEPTIDE1,1:R1-2:R4$$$V2.0", Format::HELM);
    auto linked_atomistic = toAtomistic(*linked_monomer);
    auto linked_matches = RDKit::SubstructMatch(*linked_atomistic, *tri_query);
    // 2 matches after capping group is gone (4 carbons in a chain)
    BOOST_TEST(linked_matches.size() == 2);
    BOOST_TEST(!has_atom_map_numbers(*linked_atomistic));
}

BOOST_AUTO_TEST_CASE(TestHelmImportWithCustomMonomerFromSql)
{
    // Test the SQL loading path without PDBCODE (optional field).
    // This previously crashed because _sqlite3_column_cstring returned
    // nullptr for NULL columns, causing std::string construction from
    // a null pointer.
    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.resetMonomerDefinitions();

    constexpr std::string_view custom_sql =
        "INSERT INTO monomer_definitions "
        "(SYMBOL, POLYMER_TYPE, NATURAL_ANALOG, SMILES, CORE_SMILES, "
        "NAME, MONOMER_TYPE, AUTHOR) VALUES "
        "('SqlMon', 'PEPTIDE', 'A', 'CC(C)(N[H:1])C(=O)[OH:2]', "
        "'CC(C)(N)C=O', 'SQL Monomer', 'Backbone', 'test');";

    monomer_db.loadMonomersFromSql(custom_sql);

    const std::string helm_str = "PEPTIDE1{A.[SqlMon].G}$$$$V2.0";
    auto monomer_mol = to_rdkit(helm_str, Format::HELM);
    BOOST_REQUIRE(monomer_mol != nullptr);

    auto atomistic_mol = toAtomistic(*monomer_mol);
    BOOST_REQUIRE(atomistic_mol != nullptr);
    BOOST_TEST(atomistic_mol->getNumAtoms() > 0);
}

BOOST_AUTO_TEST_CASE(TestHelmAutodetectWithUnbracketedCustomMonomer)
{
    // Multi-character monomer IDs must be bracketed per HELM spec.
    // Unbracketed should fail, bracketed should work.
    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.resetMonomerDefinitions();

    constexpr std::string_view custom_sql =
        "INSERT INTO monomer_definitions "
        "(SYMBOL, POLYMER_TYPE, NATURAL_ANALOG, SMILES, CORE_SMILES, "
        "NAME, MONOMER_TYPE, AUTHOR) VALUES "
        "('SqlMon', 'PEPTIDE', 'A', 'CC(C)(N[H:1])C(=O)[OH:2]', "
        "'CC(C)(N)C=O', 'SQL Monomer', 'Backbone', 'test');";

    monomer_db.loadMonomersFromSql(custom_sql);

    // Unbracketed multi-char monomer should fail autodetect
    const std::string unbracketed = "PEPTIDE1{A.SqlMon.G}$$$$V2.0";
    BOOST_CHECK_THROW(to_rdkit(unbracketed, Format::AUTO_DETECT),
                      std::invalid_argument);

    // Bracketed should work
    const std::string bracketed = "PEPTIDE1{A.[SqlMon].G}$$$$V2.0";
    auto monomer_mol = to_rdkit(bracketed, Format::AUTO_DETECT);
    BOOST_REQUIRE(monomer_mol != nullptr);

    auto atomistic_mol = toAtomistic(*monomer_mol);
    BOOST_REQUIRE(atomistic_mol != nullptr);
    BOOST_TEST(atomistic_mol->getNumAtoms() > 0);
}

BOOST_AUTO_TEST_CASE(TestHelmImportUnknownMonomerThrows)
{
    // Importing HELM with an unknown monomer should throw when validation
    // is enabled.
    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.resetMonomerDefinitions();

    const std::string helm_str = "PEPTIDE1{A.[UnknownMon].G}$$$$V2.0";
    BOOST_CHECK_THROW(std::ignore =
                          helm::helm_to_rdkit(helm_str, /*do_throw=*/true,
                                              /*validate=*/true),
                      std::invalid_argument);
}
