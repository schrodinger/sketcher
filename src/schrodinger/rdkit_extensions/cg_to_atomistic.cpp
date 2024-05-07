#include "schrodinger/rdkit_extensions/cg_conversions.h"

#include <queue>
#include <unordered_map>
#include <memory>
#include <fmt/format.h>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>

#include "schrodinger/rdkit_extensions/coarse_grain.h"
#include "schrodinger/rdkit_extensions/helm.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{

using AttachmentMap = std::map<std::pair<unsigned int, unsigned int>,
                               std::pair<unsigned int, unsigned int>>;

// The mapped atoms are attachment points. If used, they will be removed and
// replaced by a bond to the adjacent.
const std::unordered_map<std::string, std::string> peptide_monomer_smiles = {
    {"A", "C[C@H](N[H:1])C([OH:2])=O"},
    {"C", "[H:1]N[C@@H](CS[H:3])C([OH:2])=O"},
    {"D", "[H:1]N[C@@H](CC([OH:3])=O)C([OH:2])=O"},
    {"E", "[H:1]N[C@@H](CCC([OH:3])=O)C([OH:2])=O"},
    {"F", "[H:1]N[C@@H](Cc1ccccc1)C([OH:2])=O"},
    {"G", "[H:1]NCC([OH:2])=O"},
    {"H", "[H:1]N[C@@H](Cc1cnc[nH]1)C([OH:2])=O"},
    {"I", "CC[C@H](C)[C@H](N[H:1])C([OH:2])=O"},
    {"K", "[H:1]N[C@@H](CCCCN[H:3])C([OH:2])=O"},
    {"L", "CC(C)C[C@H](N[H:1])C([OH:2])=O"},
    {"M", "CSCC[C@H](N[H:1])C([OH:2])=O"},
    {"N", "NC(=O)C[C@H](N[H:1])C([OH:2])=O"},
    {"O", "C[C@@H]1CC=N[C@H]1C(=O)NCCCC[C@H](N[H:1])C([OH:2])=O"},
    {"P", "[H:1]N1CCC[C@H]1C([OH:2])=O"},
    {"Q", "NC(=O)CC[C@H](N[H:1])C([OH:2])=O"},
    {"R", "NC(=N)NCCC[C@H](N[H:1])C([OH:2])=O"},
    {"S", "OC[C@H](N[H:1])C([OH:2])=O"},
    {"T", "C[C@@H](O)[C@H](N[H:1])C([OH:2])=O"},
    {"U", "[SeH]C[C@H](N[H:1])C([OH:2])=O"},
    {"V", "CC(C)[C@H](N[H:1])C([OH:2])=O"},
    {"W", "[H:1]N[C@@H](Cc1c[nH]c2ccccc12)C([OH:2])=O"},
    {"Y", "Oc1ccc(C[C@H](N[H:1])C([OH:2])=O)cc1"},
    {"dA", "C[C@@H](N[H:1])C([OH:2])=O"},
    {"dC", "[H:1]N[C@H](CS[H:3])C([OH:2])=O"},
    {"dD", "[H:1]N[C@H](CC([OH:3])=O)C([OH:2])=O"},
    {"dE", "[H:1]N[C@H](CCC([OH:3])=O)C([OH:2])=O"},
    {"dF", "[H:1]N[C@H](Cc1ccccc1)C([OH:2])=O"},
    {"dH", "[H:1]N[C@H](Cc1cnc[nH]1)C([OH:2])=O"},
    {"dI", "CC[C@@H](C)[C@@H](N[H:1])C([OH:2])=O"},
    {"dK", "[H:1]N[C@H](CCCCN[H:3])C([OH:2])=O"},
    {"dL", "CC(C)C[C@@H](N[H:1])C([OH:2])=O"},
    {"dM", "CSCC[C@@H](N[H:1])C([OH:2])=O"},
    {"dN", "NC(=O)C[C@@H](N[H:1])C([OH:2])=O"},
    {"dP", "[H:1]N1CCC[C@@H]1C([OH:2])=O"},
    {"dQ", "NC(=O)CC[C@@H](N[H:1])C([OH:2])=O"},
    {"dR", "NC(=N)NCCC[C@@H](N[H:1])C([OH:2])=O"},
    {"dS", "OC[C@@H](N[H:1])C([OH:2])=O"},
    {"dT", "C[C@H](O)[C@@H](N[H:1])C([OH:2])=O"},
    {"dV", "CC(C)[C@@H](N[H:1])C([OH:2])=O"},
    {"dW", "[H:1]N[C@H](Cc1c[nH]c2ccccc12)C([OH:2])=O"},
    {"dY", "Oc1ccc(C[C@@H](N[H:1])C([OH:2])=O)cc1"},
    {"meA", "C[C@H](N(C)[H:1])C([OH:2])=O"},
    {"meC", "CN([H:1])[C@@H](CS[H:3])C([OH:2])=O"},
    {"meD", "CN([H:1])[C@@H](CC([OH:3])=O)C([OH:2])=O"},
    {"meE", "CN([H:1])[C@@H](CCC([OH:3])=O)C([OH:2])=O"},
    {"meF", "CN([H:1])[C@@H](Cc1ccccc1)C([OH:2])=O"},
    {"meG", "CN([H:1])CC([OH:2])=O"},
    {"meH", "CN([H:1])[C@@H](Cc1cnc[nH]1)C([OH:2])=O"},
    {"meI", "CC[C@H](C)[C@H](N(C)[H:1])C([OH:2])=O"},
    {"meK", "CN([H:1])[C@@H](CCCCN[H:3])C([OH:2])=O"},
    {"meL", "CC(C)C[C@H](N(C)[H:1])C([OH:2])=O"},
    {"meM", "CSCC[C@H](N(C)[H:1])C([OH:2])=O"},
    {"meN", "CN([H:1])[C@@H](CC(N)=O)C([OH:2])=O"},
    {"meQ", "CN([H:1])[C@@H](CCC(N)=O)C([OH:2])=O"},
    {"meR", "CN([H:1])[C@@H](CCCNC(N)=N)C([OH:2])=O"},
    {"meS", "CN([H:1])[C@@H](CO)C([OH:2])=O"},
    {"meT", "C[C@@H](O)[C@H](N(C)[H:1])C([OH:2])=O"},
    {"meV", "CC(C)[C@H](N(C)[H:1])C([OH:2])=O"},
    {"meW", "CN([H:1])[C@@H](Cc1c[nH]c2ccccc12)C([OH:2])=O"},
    {"meY", "CN([H:1])[C@@H](Cc1ccc(O)cc1)C([OH:2])=O"},
    {"am", "N[H:1]"},
    {"ac", "CC([OH:2])=O"}};

std::pair<unsigned int, unsigned int> get_attchpts(const std::string& linkage)
{
    // in form RX-RY, returns {X, Y}
    auto dash = linkage.find('-');
    if (dash == std::string::npos) {
        throw std::runtime_error(
            fmt::format("Invalid linkage format: {}", linkage));
    }
    return {std::stoi(linkage.substr(1, dash - 1)),
            std::stoi(linkage.substr(dash + 2))};
}

void fill_attachment_point_map(const RDKit::ROMol& new_monomer,
                               AttachmentMap& attachment_points,
                               unsigned int residue_num,
                               unsigned int old_mol_size)
{
    for (const auto& atom : new_monomer.atoms()) {
        unsigned int map_num;
        if (atom->getPropIfPresent(RDKit::common_properties::molAtomMapNumber,
                                   map_num)) {
            for (auto bnd : new_monomer.atomBonds(atom)) {
                // Should be exactly one iteration -- an attachment point can
                // only have one neighbor
                if (attachment_points.find({residue_num, map_num}) !=
                    attachment_points.end()) {
                    throw std::runtime_error(
                        fmt::format("Invalid attachment point at index {}",
                                    atom->getIdx()));
                }
                auto atom_to_bond_to =
                    old_mol_size + bnd->getOtherAtomIdx(atom->getIdx());
                auto atom_to_remove = old_mol_size + atom->getIdx();
                attachment_points[{residue_num, map_num}] =
                    std::make_pair(atom_to_bond_to, atom_to_remove);
            }
        }
    }
}

AttachmentMap add_polymer(RDKit::RWMol& atomistic_mol,
                          const RDKit::RWMol& cg_mol,
                          const std::string& polymer_id,
                          std::vector<unsigned int>& remove_atoms)
{
    // Maps residue number and attachment point number to the atom index in
    // atomistic_mol that should be attached to and the atom index of the rgroup
    // that should later be removed
    AttachmentMap attachment_point_map;

    auto chain = get_polymer(cg_mol, polymer_id);
    bool sanitize = false;

    // Add the monomers to the atomistic mol
    for (const auto monomer_idx : chain.atoms) {
        auto monomer = cg_mol.getAtomWithIdx(monomer_idx);
        auto monomer_label = monomer->getProp<std::string>(ATOM_LABEL);
        auto smiles = peptide_monomer_smiles.at(monomer_label);
        std::unique_ptr<RDKit::RWMol> new_monomer(
            RDKit::SmilesToMol(smiles, 0, sanitize));
        fill_attachment_point_map(*new_monomer, attachment_point_map,
                                  get_residue_number(monomer),
                                  atomistic_mol.getNumAtoms());
        atomistic_mol.insertMol(*new_monomer);
    }

    // Add the bonds between monomers and mark the replaced rgroups to be
    // removed
    for (const auto bond_idx : chain.bonds) {
        auto bond = cg_mol.getBondWithIdx(bond_idx);
        auto [from_rgroup, to_rgroup] =
            get_attchpts(bond->getProp<std::string>(LINKAGE));
        auto from_res = get_residue_number(bond->getBeginAtom());
        auto to_res = get_residue_number(bond->getEndAtom());
        auto [core_atom1, attachment_point1] =
            attachment_point_map.at({from_res, from_rgroup});
        auto [core_atom2, attachment_point2] =
            attachment_point_map.at({to_res, to_rgroup});
        atomistic_mol.addBond(core_atom1, core_atom2, bond->getBondType());
        remove_atoms.push_back(attachment_point1);
        remove_atoms.push_back(attachment_point2);
    }

    return attachment_point_map;
}
} // namespace

boost::shared_ptr<RDKit::RWMol> cg_to_atomistic(const RDKit::ROMol& cg_mol)
{
    boost::shared_ptr<RDKit::RWMol> atomistic_mol =
        boost::make_shared<RDKit::RWMol>();

    // Map to track Polymer ID -> attachment point map
    std::unordered_map<std::string, AttachmentMap> polymer_attachment_points;
    std::vector<unsigned int> remove_atoms;
    for (const auto& polymer_id : get_polymer_ids(cg_mol)) {
        polymer_attachment_points[polymer_id] =
            add_polymer(*atomistic_mol, cg_mol, polymer_id, remove_atoms);
    }

    // Add bonds from interpolymer connections
    for (const auto bnd : cg_mol.bonds()) {
        auto begin_atom = bnd->getBeginAtom();
        auto end_atom = bnd->getEndAtom();
        if (get_polymer_id(begin_atom) == get_polymer_id(end_atom)) {
            continue;
        }
        auto begin_res = get_residue_number(begin_atom);
        auto end_res = get_residue_number(end_atom);
        auto [from_rgroup, to_rgroup] =
            get_attchpts(bnd->getProp<std::string>(LINKAGE));
        auto [core_atom1, attachment_point1] =
            polymer_attachment_points.at(get_polymer_id(begin_atom))
                .at({begin_res, from_rgroup});
        auto [core_atom2, attachment_point2] =
            polymer_attachment_points.at(get_polymer_id(end_atom))
                .at({end_res, to_rgroup});
        atomistic_mol->addBond(core_atom1, core_atom2, bnd->getBondType());
        remove_atoms.push_back(attachment_point1);
        remove_atoms.push_back(attachment_point2);
    }

    // Remove atoms that represented attachment points
    atomistic_mol->beginBatchEdit();
    for (auto at_idx : remove_atoms) {
        atomistic_mol->removeAtom(at_idx);
    }
    atomistic_mol->commitBatchEdit();

    // Let sanitization errors bubble up for now -- it means we did something
    // wrong
    RDKit::MolOps::sanitizeMol(*atomistic_mol);

    // Remove graph hydrogens where some RGroups previously were. This will turn
    // H - NH to NH2, etc
    RDKit::MolOps::removeHs(*atomistic_mol);

    // Remove atom map numbers -- anything remaining is a capping group and the
    // map numbers are no longer meaningful
    for (auto at : atomistic_mol->atoms()) {
        at->setAtomMapNum(0);
    }

    return atomistic_mol;
}

} // namespace rdkit_extensions
} // namespace schrodinger