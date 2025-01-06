#include "schrodinger/rdkit_extensions/cg_conversions.h"

#include <queue>
#include <unordered_map>
#include <memory>

#include <boost/filesystem/path.hpp>

#include <fmt/format.h>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include "schrodinger/rdkit_extensions/cg_monomer_database.h"
#include "schrodinger/rdkit_extensions/coarse_grain.h"
#include "schrodinger/rdkit_extensions/helm.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{
namespace fs = boost::filesystem;

using AttachmentMap = std::map<std::pair<unsigned int, unsigned int>,
                               std::pair<unsigned int, unsigned int>>;

const std::unordered_map<std::string, std::string> three_character_codes({
    {"A", "ALA"}, // Alanine
    {"R", "ARG"}, // Arginine
    {"N", "ASN"}, // Asparagine
    {"D", "ASP"}, // Aspartic
    {"C", "CYS"}, // Cysteine
    {"Q", "GLN"}, // Glutamine
    {"E", "GLU"}, // Glutamic
    {"G", "GLY"}, // Glycine
    {"H", "HIS"}, // Histidine
    {"I", "ILE"}, // Isoleucine
    {"L", "LEU"}, // Leucine
    {"K", "LYS"}, // Lysine
    {"M", "MET"}, // Methionine
    {"F", "PHE"}, // Phenylalanine
    {"P", "PRO"}, // Proline
    {"S", "SER"}, // Serine
    {"T", "THR"}, // Threonine
    {"W", "TRP"}, // Tryptophan
    {"Y", "TYR"}, // Tyrosine
    {"V", "VAL"}, // Valine
});

static const std::string ATOM_PDB_NAME_PROP{"pdbName"};

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

void set_pdb_info(RDKit::RWMol& new_monomer, const std::string& monomer_label,
                  unsigned int residue_number, char chain_id,
                  ChainType chain_type)
{
    std::string residue_name =
        (chain_type == ChainType::PEPTIDE) ? "UNK" : "UNL";
    if (three_character_codes.find(monomer_label) !=
        three_character_codes.end()) {
        residue_name = three_character_codes.at(monomer_label);
    }

    // Set PDB residue info on the atoms
    for (auto atom : new_monomer.atoms()) {
        auto* res_info = new RDKit::AtomPDBResidueInfo();

        std::string chain_id_str(1, chain_id);
        res_info->setChainId(chain_id_str);
        res_info->setResidueNumber(residue_number);
        res_info->setResidueName(residue_name);
        // to be consistent with the rdkit adapter and RDKit's own PDB writer
        res_info->setInsertionCode(" ");

        std::string pdb_name;
        if (atom->getPropIfPresent(ATOM_PDB_NAME_PROP, pdb_name)) {
            res_info->setName(pdb_name);
            // Don't keep property on output molecule
            atom->clearProp(ATOM_PDB_NAME_PROP);
        }

        atom->setMonomerInfo(res_info);
    }
}

ChainType get_chain_type(std::string_view polymer_id)
{
    if (polymer_id.find("PEPTIDE") == 0) {
        return ChainType::PEPTIDE;
    } else if (polymer_id.find("RNA") == 0) {
        // HELM labels both DNA and RNA as RNA
        return ChainType::RNA;
    } else if (polymer_id.find("CHEM") == 0) {
        return ChainType::CHEM;
    } else {
        throw std::out_of_range(fmt::format(
            "Invalid polymer id: {}. Must be one of PEPTIDE, RNA, CHEM",
            polymer_id));
    }
}

AttachmentMap add_polymer(RDKit::RWMol& atomistic_mol,
                          const RDKit::RWMol& cg_mol,
                          const std::string& polymer_id,
                          std::vector<unsigned int>& remove_atoms,
                          char chain_id)
{
    // Maps residue number and attachment point number to the atom index in
    // atomistic_mol that should be attached to and the atom index of the rgroup
    // that should later be removed
    AttachmentMap attachment_point_map;

    auto chain = get_polymer(cg_mol, polymer_id);
    auto chain_type = get_chain_type(polymer_id);
    bool sanitize = false;

    cg_monomer_database db(get_cg_monomer_db_path());

    // Add the monomers to the atomistic mol
    for (const auto monomer_idx : chain.atoms) {
        auto monomer = cg_mol.getAtomWithIdx(monomer_idx);
        auto monomer_label = monomer->getProp<std::string>(ATOM_LABEL);

        std::string smiles;
        if (monomer->getProp<bool>(SMILES_MONOMER)) {
            smiles = monomer_label;
        } else {
            auto monomer_smiles =
                db.get_monomer_smiles(monomer_label, chain_type);
            if (!monomer_smiles) {
                throw std::out_of_range(fmt::format(
                    "Peptide Monomer {} not found in CG Monomer database",
                    monomer_label));
            }
            smiles = *monomer_smiles;
        }

        std::unique_ptr<RDKit::RWMol> new_monomer(
            RDKit::SmilesToMol(smiles, 0, sanitize));

        if (monomer->getProp<bool>(SMILES_MONOMER)) {
            // SMILES monomers may be in rgroup form like
            // *N[C@H](C(=O)O)S* |$_R1;;;;;;;_R3$| or use atom map numbers like
            // [*:1]N[C@H](C(=O)O)S[*:3]. Translate the RGroup to atom map
            // numbers
            for (auto atom : new_monomer->atoms()) {
                std::string rgroup_label;
                if (atom->getPropIfPresent(RDKit::common_properties::atomLabel,
                                           rgroup_label) &&
                    rgroup_label.find("_R") == 0) {
                    auto rgroup_num = std::stoi(rgroup_label.substr(2));
                    atom->setAtomMapNum(rgroup_num);
                    atom->clearProp(RDKit::common_properties::atomLabel);
                }
            }
        }

        auto residue_number = get_residue_number(monomer);
        fill_attachment_point_map(*new_monomer, attachment_point_map,
                                  residue_number, atomistic_mol.getNumAtoms());
        set_pdb_info(*new_monomer, monomer_label, residue_number, chain_id,
                     chain_type);

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

        if (attachment_point_map.find({from_res, from_rgroup}) ==
                attachment_point_map.end() ||
            attachment_point_map.find({to_res, to_rgroup}) ==
                attachment_point_map.end()) {
            // One of these attachment points is not present
            throw std::runtime_error(fmt::format(
                "Invalid linkage {} between monomers {} and {}",
                bond->getProp<std::string>(LINKAGE), from_res, to_res));
        }

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
    auto atomistic_mol = boost::make_shared<RDKit::RWMol>();

    // Map to track Polymer ID -> attachment point map
    std::unordered_map<std::string, AttachmentMap> polymer_attachment_points;
    std::vector<unsigned int> remove_atoms;
    char chain_id = 'A';
    for (const auto& polymer_id : get_polymer_ids(cg_mol)) {
        polymer_attachment_points[polymer_id] = add_polymer(
            *atomistic_mol, cg_mol, polymer_id, remove_atoms, chain_id);
        ++chain_id;
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

        const auto& begin_attachment_points =
            polymer_attachment_points.at(get_polymer_id(begin_atom));
        const auto& end_attachment_points =
            polymer_attachment_points.at(get_polymer_id(end_atom));

        if (begin_attachment_points.find({begin_res, from_rgroup}) ==
                begin_attachment_points.end() ||
            end_attachment_points.find({end_res, to_rgroup}) ==
                end_attachment_points.end()) {
            // One of these attachment points is not present
            throw std::runtime_error(
                fmt::format("Invalid linkage {} between monomers {} and {}",
                            bnd->getProp<std::string>(LINKAGE),
                            begin_atom->getIdx(), end_atom->getIdx()));
        }

        auto [core_atom1, attachment_point1] =
            begin_attachment_points.at({begin_res, from_rgroup});
        auto [core_atom2, attachment_point2] =
            end_attachment_points.at({end_res, to_rgroup});
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
