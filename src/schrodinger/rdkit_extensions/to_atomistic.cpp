#include "schrodinger/rdkit_extensions/atomistic_conversions.h"

#include <algorithm>
#include <cctype>
#include <memory>
#include <optional>
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/algorithm/string/trim.hpp>
#include <boost/json.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/system/system_error.hpp>

#include <fmt/format.h>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>
#include <rdkit/GraphMol/Substruct/SubstructMatch.h>
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/sup_utils.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{
namespace fs = boost::filesystem;

const std::string ATOM_PROP_GRAPH_H("_GRAPH_H_TO_REMOVE");
const std::string HYDROGEN_BOND_LINKAGE{"pair-pair"};

struct AttachmentInfo {
    unsigned int core_atom;
    std::vector<unsigned int> leaving_group_atoms;
};

using AttachmentMap =
    std::map<std::pair<unsigned int, unsigned int>, AttachmentInfo>;
using AtomIndexSet = std::unordered_set<unsigned int>;

static const std::string ATOM_PDB_NAME_PROP{"pdbName"};

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

static const std::map<unsigned int, std::string> BIOVIA_ATTCHPT_MAP = {
    {1, "Al"}, // Backbone attachment point
    {2, "Br"}, // Backbone attachment point
    {3, "Cx"}, // Sidechain attachment point
};

struct MonomerKey {
    std::string chain_id;
    unsigned int residue_num;

    bool operator<(const MonomerKey& other) const
    {
        return std::tie(chain_id, residue_num) <
               std::tie(other.chain_id, other.residue_num);
    }
};

const std::map<std::string, std::vector<std::pair<std::string, std::string>>>
    BASE_PAIRWISE_HBONDS = {
        // Adenine - Thymine (2 bonds)
        {"AT", {{"N1", "N3"}, {"N6", "O4"}}},
        {"TA", {{"N3", "N1"}, {"O4", "N6"}}},

        // Adenine - Uracil (2 bonds)
        {"AU", {{"N1", "N3"}, {"N6", "O4"}}},
        {"UA", {{"N3", "N1"}, {"O4", "N6"}}},

        // Guanine - Cytosine (3 bonds)
        {"GC", {{"O6", "N4"}, {"N1", "N3"}, {"N2", "O2"}}},
        {"CG", {{"N4", "O6"}, {"N3", "N1"}, {"O2", "N2"}}}};

std::pair<unsigned int, unsigned int> getAttchpts(const std::string& linkage)
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

std::optional<unsigned int>
findMappedAttachmentAtom(const RDKit::ROMol& monomer, unsigned int rgroup_num)
{
    for (const auto* atom : monomer.atoms()) {
        unsigned int map_num = 0;
        if (atom->getPropIfPresent(RDKit::common_properties::molAtomMapNumber,
                                   map_num) &&
            map_num == rgroup_num) {
            return atom->getIdx();
        }
    }
    return std::nullopt;
}

AtomIndexSet matchCoreAtoms(const RDKit::ROMol& new_monomer,
                            const std::string& core_smiles)
{
    std::unique_ptr<RDKit::RWMol> core_mol(
        RDKit::SmilesToMol(core_smiles, 0, false));
    if (!core_mol) {
        throw std::runtime_error(fmt::format(
            "Invalid CORE_SMILES for R-group matching: {}", core_smiles));
    }

    // Attachment-point dummies describe bonds leaving the core. Keep the core
    // atoms as the query and remove the dummies before matching.
    std::vector<unsigned int> dummy_indices;
    for (const auto* atom : core_mol->atoms()) {
        if (atom->getAtomicNum() == 0) {
            dummy_indices.push_back(atom->getIdx());
        }
    }
    std::sort(dummy_indices.rbegin(), dummy_indices.rend());
    for (auto atom_idx : dummy_indices) {
        core_mol->removeAtom(atom_idx);
    }
    if (core_mol->getNumAtoms() == 0) {
        throw std::runtime_error(
            fmt::format("CORE_SMILES contains no core atoms: {}", core_smiles));
    }

    const auto core_matches = RDKit::SubstructMatch(new_monomer, *core_mol);
    if (core_matches.empty()) {
        throw std::runtime_error(fmt::format(
            "Could not match CORE_SMILES in monomer: {}", core_smiles));
    }

    // Multiple query mappings can describe the same set of core atoms (for
    // example, due to symmetry). They define the same core/non-core split.
    AtomIndexSet core_atoms;
    for (const auto& [query_idx, target_idx] : core_matches.front()) {
        (void)query_idx;
        core_atoms.insert(static_cast<unsigned int>(target_idx));
    }
    for (const auto& core_match : core_matches) {
        AtomIndexSet matched_core_atoms;
        for (const auto& [query_idx, target_idx] : core_match) {
            (void)query_idx;
            matched_core_atoms.insert(static_cast<unsigned int>(target_idx));
        }
        if (matched_core_atoms != core_atoms) {
            throw std::runtime_error(
                fmt::format("CORE_SMILES matched multiple core atom sets: {}",
                            core_smiles));
        }
    }
    return core_atoms;
}

boost::json::value parseRGroups(const std::string& rgroups)
{
    boost::json::value json_rgroups;
    try {
        json_rgroups = boost::json::parse(rgroups);
    } catch (const boost::system::system_error& e) {
        throw std::runtime_error(
            fmt::format("Invalid RGROUPS JSON: {}", e.what()));
    }
    return json_rgroups;
}

unsigned int parseRGroup(boost::json::string_view rgroup_name,
                         const boost::json::value& rgroup_value)
{
    if (!rgroup_value.is_string()) {
        throw std::runtime_error(fmt::format(
            "RGROUPS entry {} must contain a SMILES string", rgroup_name));
    }

    auto rgroup_mol = std::unique_ptr<RDKit::RWMol>(
        RDKit::SmilesToMol(rgroup_value.as_string().c_str(), 0, false));
    if (!rgroup_mol) {
        throw std::runtime_error(
            fmt::format("Invalid SMILES for R-group {}", rgroup_name));
    }

    RDKit::Atom* rgroup_dummy_atom = nullptr;
    for (auto* atom : rgroup_mol->atoms()) {
        if (atom->getAtomicNum() == 0) {
            if (rgroup_dummy_atom != nullptr) {
                throw std::runtime_error(fmt::format(
                    "R-group {} contains multiple dummy atoms", rgroup_name));
            }
            rgroup_dummy_atom = atom;
        }
    }
    if (rgroup_dummy_atom == nullptr) {
        throw std::runtime_error(fmt::format(
            "R-group {} does not contain a dummy atom", rgroup_name));
    }

    // The rgroup num matches the atom map num in the mol
    const auto rgroup_num = static_cast<unsigned int>(
        std::stoul(std::string(rgroup_name.substr(1))));
    rgroup_mol->removeAtom(rgroup_dummy_atom->getIdx());
    if (rgroup_mol->getNumAtoms() == 0) {
        throw std::runtime_error(
            fmt::format("R-group {} is empty", rgroup_name));
    }
    return rgroup_num;
}

AtomIndexSet findLeavingGroupAtoms(const RDKit::ROMol& new_monomer,
                                   const AtomIndexSet& core_atoms,
                                   unsigned int mapped_attachment_atom)
{
    // Traverse the non-core component starting from the R-group's mapped
    // atom. The core atom set acts as a boundary, so similar groups
    // elsewhere in the monomer cannot be mistaken for this leaving group.
    AtomIndexSet lg_atoms;
    std::queue<unsigned int> atoms_to_visit;
    lg_atoms.insert(mapped_attachment_atom);
    atoms_to_visit.push(mapped_attachment_atom);
    while (!atoms_to_visit.empty()) {
        const auto atom_idx = atoms_to_visit.front();
        atoms_to_visit.pop();
        const auto* atom = new_monomer.getAtomWithIdx(atom_idx);
        for (const auto* neighbor : new_monomer.atomNeighbors(atom)) {
            const auto neighbor_idx = neighbor->getIdx();
            if (core_atoms.count(neighbor_idx) == 0 &&
                lg_atoms.insert(neighbor_idx).second) {
                atoms_to_visit.push(neighbor_idx);
            }
        }
    }
    return lg_atoms;
}

unsigned int findAttachedCoreAtom(const RDKit::ROMol& new_monomer,
                                  const AtomIndexSet& lg_atoms,
                                  boost::json::string_view rgroup_name)
{
    AtomIndexSet attach_core_atoms;
    for (auto group_atom_idx : lg_atoms) {
        auto* atom = new_monomer.getAtomWithIdx(group_atom_idx);
        for (auto* neighbor : new_monomer.atomNeighbors(atom)) {
            if (!lg_atoms.contains(neighbor->getIdx())) {
                attach_core_atoms.insert(neighbor->getIdx());
            }
        }
    }
    // The LG should be attached to 1 core atom
    if (attach_core_atoms.size() != 1) {
        throw std::runtime_error(fmt::format(
            "Could not identify the attached core atom for R-group {}",
            rgroup_name));
    }
    return *attach_core_atoms.begin();
}

void fillAttachmentPointMapFromRGroups(const RDKit::ROMol& new_monomer,
                                       const std::string& rgroups,
                                       const std::string& core_smiles,
                                       AttachmentMap& attachment_points,
                                       unsigned int residue_num,
                                       unsigned int old_mol_size)
{
    // Match CORE_SMILES first. The atoms outside this core, reachable from
    // each R-group's mapped atom, make up that R-group's leaving group.
    const auto core_atoms = matchCoreAtoms(new_monomer, core_smiles);
    const auto json_rgroups = parseRGroups(rgroups);
    for (const auto& [rgroup_name, rgroup_value] : json_rgroups.as_object()) {
        const auto rgroup_num = parseRGroup(rgroup_name, rgroup_value);
        const auto mapped_attachment_atom =
            findMappedAttachmentAtom(new_monomer, rgroup_num);
        if (!mapped_attachment_atom) {
            throw std::runtime_error(
                fmt::format("Could not find mapped atom for R-group {} in "
                            "monomer SMILES",
                            rgroup_name));
        }
        if (core_atoms.count(*mapped_attachment_atom) != 0) {
            throw std::runtime_error(
                fmt::format("Mapped atom for R-group {} is part of CORE_SMILES",
                            rgroup_name));
        }

        const auto lg_atoms = findLeavingGroupAtoms(new_monomer, core_atoms,
                                                    *mapped_attachment_atom);
        const auto attached_core_atom =
            findAttachedCoreAtom(new_monomer, lg_atoms, rgroup_name);

        // Offset local monomer atom indicies for the overall polymer
        AttachmentInfo info{old_mol_size + attached_core_atom, {}};
        info.leaving_group_atoms.reserve(lg_atoms.size());
        for (auto atom_idx : lg_atoms) {
            info.leaving_group_atoms.push_back(old_mol_size + atom_idx);
        }
        attachment_points[{residue_num, rgroup_num}] = std::move(info);
    }
}

void fillAttachmentPointMapFromSmiles(const RDKit::ROMol& new_monomer,
                                      AttachmentMap& attachment_points,
                                      unsigned int residue_num,
                                      unsigned int old_mol_size)
{
    // Fallback method to determine attachment point map using SMILES
    // Only 1 heavy atom LGs are supported, any LG with more than 1
    // heavy atom neighbor will throw an error.
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
                        fmt::format("Multiple heavy atom LG for attachment "
                                    "point on atom {}, use RGROUPS + "
                                    "CORE_SMILES instead.",
                                    atom->getIdx()));
                }
                auto atom_to_bond_to =
                    old_mol_size + bnd->getOtherAtomIdx(atom->getIdx());
                auto atom_to_remove = old_mol_size + atom->getIdx();
                attachment_points[{residue_num, map_num}] = {atom_to_bond_to,
                                                             {atom_to_remove}};
            }
        }
    }
}

void fillAttachmentPointMap(const RDKit::ROMol& new_monomer,
                            const MonomerInfo* monomer_info,
                            AttachmentMap& attachment_points,
                            unsigned int residue_num, unsigned int old_mol_size)
{
    if (monomer_info && monomer_info->rgroups) {
        auto rgroups = boost::trim_copy(*monomer_info->rgroups);
        // The default database represents an empty RGROUPS field as the
        // serialized empty JSON string, rather than as an empty field.
        if (!rgroups.empty() && rgroups != "\"\"") {
            if (!monomer_info->core_smiles ||
                monomer_info->core_smiles->empty()) {
                throw std::runtime_error(
                    "CORE_SMILES is required when RGROUPS are defined");
            }
            fillAttachmentPointMapFromRGroups(
                new_monomer, rgroups, *monomer_info->core_smiles,
                attachment_points, residue_num, old_mol_size);
            return;
        }
    }

    fillAttachmentPointMapFromSmiles(new_monomer, attachment_points,
                                     residue_num, old_mol_size);
}

void setResidueInfo(RDKit::RWMol& new_monomer, const std::string& monomer_label,
                    unsigned int residue_number, char chain_id,
                    ChainType chain_type, unsigned int current_residue,
                    const MonomerDatabase& db)
{
    std::string residue_name =
        (chain_type == ChainType::PEPTIDE) ? "UNK" : "UNL";

    auto pdb_code = db.getPdbCode(monomer_label, chain_type);
    if (pdb_code) {
        residue_name = *pdb_code;
    } else if (three_character_codes.find(monomer_label) !=
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

std::vector<std::pair<unsigned int, unsigned int>> getHBondAtomIndices(
    const std::vector<std::pair<std::string, std::string>>& hbonds,
    const std::vector<RDKit::Atom*>& begin_base_atoms,
    const std::vector<RDKit::Atom*>& end_base_atoms)
{
    std::vector<std::pair<unsigned int, unsigned int>> hbond_pair_indices;
    for (auto& hbond : hbonds) {
        int first_idx = -1;
        int second_idx = -1;

        for (auto* atom : begin_base_atoms) {
            auto* pdb_info = dynamic_cast<RDKit::AtomPDBResidueInfo*>(
                atom->getMonomerInfo());
            std::string atom_name = pdb_info->getName();
            boost::trim(atom_name);
            if (atom_name == hbond.first) {
                first_idx = atom->getIdx();
                break;
            }
        }
        for (auto* atom : end_base_atoms) {
            auto* pdb_info = dynamic_cast<RDKit::AtomPDBResidueInfo*>(
                atom->getMonomerInfo());
            std::string atom_name = pdb_info->getName();
            boost::trim(atom_name);
            if (atom_name == hbond.second) {
                second_idx = atom->getIdx();
                break;
            }
        }

        // Validate that both atom indices were found
        if (first_idx == -1 || second_idx == -1) {
            throw std::runtime_error(
                fmt::format("Could not find hydrogen bond atoms: {} - {}",
                            hbond.first, hbond.second));
        }
        hbond_pair_indices.push_back({static_cast<unsigned int>(first_idx),
                                      static_cast<unsigned int>(second_idx)});
    }
    return hbond_pair_indices;
}

AttachmentMap addPolymer(RDKit::RWMol& atomistic_mol,
                         const RDKit::RWMol& monomer_mol,
                         const std::string& polymer_id,
                         std::vector<unsigned int>& remove_atoms, char chain_id,
                         unsigned int& total_residue_count)
{
    // Maps residue number and attachment point number to the atom index in
    // atomistic_mol that should be attached to and the atom index of the rgroup
    // that should later be removed
    AttachmentMap attachment_point_map;
    auto chain = get_polymer(monomer_mol, polymer_id);
    auto chain_type = getChainType(polymer_id);
    bool sanitize = false;

    auto& db = MonomerDatabase::instance();

    // Add the monomers to the atomistic mol
    for (const auto monomer_idx : chain.atoms) {
        auto monomer = monomer_mol.getAtomWithIdx(monomer_idx);
        auto monomer_label = monomer->getProp<std::string>(ATOM_LABEL);
        auto monomer_info = db.getMonomerInfo(monomer_label, chain_type);

        std::string smiles;
        if (monomer->getProp<bool>(SMILES_MONOMER)) {
            smiles = monomer_label;
        } else {
            auto monomer_smiles =
                db.getMonomerSmiles(monomer_label, chain_type);
            if (!monomer_smiles) {
                throw std::out_of_range(fmt::format(
                    "Peptide Monomer {} not found in Monomer database",
                    monomer_label));
            }
            smiles = *monomer_smiles;
        }

        std::unique_ptr<RDKit::RWMol> new_monomer(
            RDKit::SmilesToMol(smiles, 0, sanitize));

        if (!new_monomer) {
            // FIXME: I think this is an issue with the HELM parser, see
            // SHARED-11457
            new_monomer.reset(
                RDKit::SmilesToMol("[" + smiles + "]", 0, sanitize));
        }

        if (!new_monomer) {
            throw std::runtime_error(
                fmt::format("Monomer has invalid SMILES: {}", smiles));
        }

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

        // Mark graph hydrogens with atom map numbers for later removal
        for (auto atom : new_monomer->atoms()) {
            if (atom->getAtomicNum() == 1) {
                unsigned int map_num;
                if (atom->getPropIfPresent(
                        RDKit::common_properties::molAtomMapNumber, map_num)) {
                    atom->setProp(ATOM_PROP_GRAPH_H, true);
                }
            }
        }

        auto residue_number = get_residue_number(monomer);
        fillAttachmentPointMap(
            *new_monomer,
            monomer->getProp<bool>(SMILES_MONOMER) ? nullptr : &monomer_info,
            attachment_point_map, residue_number, atomistic_mol.getNumAtoms());
        setResidueInfo(*new_monomer, monomer_label, residue_number, chain_id,
                       chain_type, total_residue_count, db);
        ++total_residue_count;
        atomistic_mol.insertMol(*new_monomer);
    }

    std::map<unsigned int, std::vector<RDKit::Atom*>> atoms_by_monomer;
    for (auto atom : atomistic_mol.atoms()) {
        auto chain = get_polymer_id(atom);
        // Only get the residues for this polymer
        if (chain == std::string(1, chain_id)) {
            auto res_num = get_residue_number(atom);
            atoms_by_monomer[res_num].push_back(atom);
        }
    }
    // Add the bonds between monomers and mark the replaced rgroups to be
    // removed
    for (const auto bond_idx : chain.bonds) {
        auto bond = monomer_mol.getBondWithIdx(bond_idx);
        auto linkage = bond->getProp<std::string>(LINKAGE);
        auto begin_atom = bond->getBeginAtom();
        auto end_atom = bond->getEndAtom();
        auto from_res = get_residue_number(begin_atom);
        auto to_res = get_residue_number(end_atom);
        if (linkage == HYDROGEN_BOND_LINKAGE) {
            auto from_base = begin_atom->getProp<std::string>(ATOM_LABEL);
            auto begin_base_atoms = atoms_by_monomer[from_res];

            auto end_base = end_atom->getProp<std::string>(ATOM_LABEL);
            auto end_base_atoms = atoms_by_monomer[to_res];

            auto pair = from_base + end_base;
            auto it = BASE_PAIRWISE_HBONDS.find(pair);
            if (it == BASE_PAIRWISE_HBONDS.end()) {
                throw std::runtime_error(
                    fmt::format("Hydrogen Bonds between the pair of bases "
                                "{}-{} is not supported.",
                                from_base, end_base));
            }
            auto hbond_pair_indices = getHBondAtomIndices(
                it->second, begin_base_atoms, end_base_atoms);
            for (auto& hbond : hbond_pair_indices) {
                // The bond-type is ZERO
                atomistic_mol.addBond(hbond.first, hbond.second,
                                      bond->getBondType());
            }
            continue;
        }
        // Add custom bonds to linkages to add in the atomisic rep
        std::vector<std::string> linkages = {linkage};
        std::string custom_bond;
        if (bond->getPropIfPresent<std::string>(CUSTOM_BOND, custom_bond)) {
            if (linkage != custom_bond) {
                linkages.push_back(custom_bond);
            }
        }
        for (auto& link : linkages) {
            auto [from_rgroup, to_rgroup] = getAttchpts(link);
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

            auto [core_aid1, attachment_info1] =
                attachment_point_map.at({from_res, from_rgroup});
            auto [core_aid2, attachment_info2] =
                attachment_point_map.at({to_res, to_rgroup});

            auto bond_type = bond->getBondType();
            if (bond_type == RDKit::Bond::DATIVE) {
                // Only relevant at the monomer mol level, this is just
                // a single bond at the atomistic level
                bond_type = RDKit::Bond::SINGLE;
            }
            atomistic_mol.addBond(core_aid1, core_aid2, bond_type);
            remove_atoms.insert(remove_atoms.end(), attachment_info1.begin(),
                                attachment_info1.end());
            remove_atoms.insert(remove_atoms.end(), attachment_info2.begin(),
                                attachment_info2.end());
        }
    }

    return attachment_point_map;
}
} // namespace

boost::shared_ptr<RDKit::RWMol> toAtomistic(const RDKit::ROMol& monomer_mol)
{
    auto atomistic_mol = boost::make_shared<RDKit::RWMol>();

    // Map to track Polymer ID -> attachment point map
    std::unordered_map<std::string, AttachmentMap> polymer_attachment_points;
    std::vector<unsigned int> remove_atoms;
    unsigned int total_residue_count = 1; // 1-based index to label SUP groups
    char chain_id = 'A';
    std::map<std::string, std::string> polymer_to_chain;
    for (const auto& polymer_id : get_polymer_ids(monomer_mol)) {
        polymer_attachment_points[polymer_id] =
            addPolymer(*atomistic_mol, monomer_mol, polymer_id, remove_atoms,
                       chain_id, total_residue_count);
        polymer_to_chain[polymer_id] = chain_id;
        ++chain_id;
    }

    std::map<MonomerKey, std::vector<RDKit::Atom*>> atoms_by_monomer;
    for (auto atom : atomistic_mol->atoms()) {
        // In atomistic mol polymer id is chain
        auto chain = get_polymer_id(atom);
        auto res_num = get_residue_number(atom);
        MonomerKey mono_key{chain, res_num};
        atoms_by_monomer[mono_key].push_back(atom);
    }

    // Add bonds from interpolymer connections
    for (const auto bnd : monomer_mol.bonds()) {
        auto begin_atom = bnd->getBeginAtom();
        auto end_atom = bnd->getEndAtom();
        if (get_polymer_id(begin_atom) == get_polymer_id(end_atom)) {
            continue;
        }
        auto begin_res = get_residue_number(begin_atom);
        auto end_res = get_residue_number(end_atom);
        auto linkage = bnd->getProp<std::string>(LINKAGE);
        if (linkage == HYDROGEN_BOND_LINKAGE) {
            auto from_base = begin_atom->getProp<std::string>(ATOM_LABEL);
            MonomerKey begin_base_key{
                polymer_to_chain[get_polymer_id(begin_atom)], begin_res};
            auto begin_base_atoms = atoms_by_monomer[begin_base_key];

            auto end_base = end_atom->getProp<std::string>(ATOM_LABEL);
            MonomerKey end_base_key{polymer_to_chain[get_polymer_id(end_atom)],
                                    end_res};
            auto end_base_atoms = atoms_by_monomer[end_base_key];

            auto pair = from_base + end_base;
            auto it = BASE_PAIRWISE_HBONDS.find(pair);
            if (it == BASE_PAIRWISE_HBONDS.end()) {
                throw std::runtime_error(
                    fmt::format("Hydrogen Bonds between the pair of bases "
                                "{}-{} is not supported.",
                                from_base, end_base));
            }
            auto hbond_pair_indices = getHBondAtomIndices(
                it->second, begin_base_atoms, end_base_atoms);
            for (auto& hbond : hbond_pair_indices) {
                // The bond-type is ZERO
                atomistic_mol->addBond(hbond.first, hbond.second,
                                       bnd->getBondType());
            }
            continue;
        }
        auto [from_rgroup, to_rgroup] = getAttchpts(linkage);

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

        auto [core_atom1, attachment_info1] =
            begin_attachment_points.at({begin_res, from_rgroup});
        auto [core_atom2, attachment_info2] =
            end_attachment_points.at({end_res, to_rgroup});
        auto bond_type = bnd->getBondType();
        if (bond_type == RDKit::Bond::DATIVE) {
            // Only relevant at the monomer mol level, this is just
            // a single bond at the atomistic level
            bond_type = RDKit::Bond::SINGLE;
        }
        atomistic_mol->addBond(core_atom1, core_atom2, bond_type);
        remove_atoms.insert(remove_atoms.end(), attachment_info1.begin(),
                            attachment_info1.end());
        remove_atoms.insert(remove_atoms.end(), attachment_info2.begin(),
                            attachment_info2.end());
    }

    // Remove atoms that represented attachment points/dummy atoms
    // and graph hydrogens marked earlier. The latter is faster than
    // calling removeHs explicitly which checks all hydrogens.
    atomistic_mol->beginBatchEdit();
    for (auto at_idx : remove_atoms) {
        auto* at = atomistic_mol->getAtomWithIdx(at_idx);
        atomistic_mol->removeAtom(at_idx);

        // We already removed the atom because it was used as an attachment
        // point, so we don't want to remove it again below.
        at->clearProp(ATOM_PROP_GRAPH_H);
    }
    for (auto at : atomistic_mol->atoms()) {
        bool is_graph_h = false;
        if (at->getPropIfPresent(ATOM_PROP_GRAPH_H, is_graph_h) && is_graph_h) {
            // Inc. neighbor's H count to prevent pyrrole kekulization issues.
            for (auto nbr : atomistic_mol->atomNeighbors(at)) {
                nbr->setNumExplicitHs(nbr->getNumExplicitHs() + 1);
            }
            atomistic_mol->removeAtom(at->getIdx());
        } else if (at->getAtomicNum() == 0) {
            atomistic_mol->removeAtom(at->getIdx());
        }
    }
    atomistic_mol->commitBatchEdit();

    // Let sanitization errors bubble up for now -- it means we did something
    // wrong
    RDKit::MolOps::sanitizeMol(*atomistic_mol);

    // Remove atom map numbers -- anything remaining is a capping group and the
    // map numbers are no longer meaningful
    for (auto at : atomistic_mol->atoms()) {
        at->setAtomMapNum(0);
    }

    // Add SUP substance groups to annotate residues for SDFs
    createSupGroupsFromResidueInfo(*atomistic_mol);

    CopyMolProperties(monomer_mol, *atomistic_mol);
    return atomistic_mol;
}

} // namespace rdkit_extensions
} // namespace schrodinger
