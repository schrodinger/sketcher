#include "schrodinger/rdkit_extensions/sup_utils.h"

#include <algorithm>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include "schrodinger/rdkit_extensions/monomer_database.h"

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{

const std::string MONOMER_IDX1{"monomerIndex1"};

bool hasSupSgroups(const RDKit::ROMol& mol)
{
    // Check for the presence of SUP substance groups
    for (const auto& sgroup : RDKit::getSubstanceGroups(mol)) {
        std::string type;
        if (sgroup.getPropIfPresent(SGROUP_PROP_TYPE, type) && type == "SUP") {
            return true;
        }
    }
    return false;
}

void supGroupsToPdbInfo(RDKit::ROMol& mol)
{
    // If SUP groups aren't labeled in any order, set residue number
    const auto& sgroups = RDKit::getSubstanceGroups(mol);
    const auto it = std::find_if(sgroups.begin(), sgroups.end(),
                                 [](const RDKit::SubstanceGroup& sgroup) {
                                     std::string type;
                                     return sgroup.getPropIfPresent(
                                                SGROUP_PROP_TYPE, type) &&
                                            type == "SUP";
                                 });
    if (it == sgroups.end()) {
        throw std::runtime_error("No SUP sgroups found");
    }

    if (!sgroups.empty() && !(*it).hasProp(SGROUP_PROP_RESNUM)) {
        unsigned int res_num = 1;
        for (const auto& sgroup : sgroups) {
            sgroup.setProp<unsigned int>(SGROUP_PROP_RESNUM, res_num);
            ++res_num;
        }
    }

    // SUP sgroups, or abbreviations, are sometimes set to indicate monomers
    // on an atomistic molecule. They do not indicate any chain information,
    // so for now we will place them all in a single chain.
    std::string chain_id = "A";

    // Convert SUP sgroups to PDB residue info
    for (const auto& sgroup : sgroups) {
        std::string type;
        if (sgroup.getPropIfPresent(SGROUP_PROP_TYPE, type) && type == "SUP") {
            for (const auto& atom_idx : sgroup.getAtoms()) {
                auto atom = mol.getAtomWithIdx(atom_idx);
                std::string symbol;
                std::string sup_class;
                unsigned int res_num = 0;
                if (!sgroup.getPropIfPresent(SGROUP_PROP_LABEL, symbol) ||
                    !sgroup.getPropIfPresent(SGROUP_PROP_RESNUM, res_num) ||
                    !sgroup.getPropIfPresent(SGROUP_PROP_CLASS, sup_class)) {
                    throw std::runtime_error(
                        "SUP sgroup missing label or class");
                }
                auto res_info = new RDKit::AtomPDBResidueInfo();
                res_info->setResidueName(symbol);
                res_info->setResidueNumber(res_num);
                res_info->setChainId(chain_id);
                atom->setMonomerInfo(res_info);
                atom->setProp<std::string>(SGROUP_PROP_CLASS, sup_class);
            }
        }
    }
}

void removeHydrogenSupGroups(RDKit::ROMol& mol)
{
    auto& sgroups = RDKit::getSubstanceGroups(mol);
    // Track hydrogens that are "apart of the chain" (backbone)
    for (const auto& sg : sgroups) {
        for (auto& attch : sg.getAttachPoints()) {
            if ((attch.id == "Al" || attch.id == "Br") &&
                mol.getAtomWithIdx(attch.lvIdx)->getAtomicNum() == 1) {
                mol.getAtomWithIdx(attch.lvIdx)->setProp(CHAIN_HYDROGEN, true);
            }
        }
    }
    // Remove any hydrogen sgroups from the molecule that isn't in the
    // "backbone"
    std::vector<RDKit::SubstanceGroup> newsgs;
    for (auto& sgroup : RDKit::getSubstanceGroups(mol)) {
        std::string label;
        auto at = mol.getAtomWithIdx(sgroup.getAtoms()[0]);
        if (sgroup.getPropIfPresent(SGROUP_PROP_LABEL, label) && label == "H" &&
            sgroup.getAtoms().size() == 1 && !at->hasProp(CHAIN_HYDROGEN)) {
            continue;
        } else {
            newsgs.push_back(std::move(sgroup));
        }
    }
    sgroups = std::move(newsgs);
}

void addUnmappedAtomsToResidues(RDKit::ROMol& mol)
{
    // Add unmapped atoms to the residue they belong to
    // make RWMol
    RDKit::RWMol rwmol(mol);
    std::vector<unsigned int> to_remove;
    int res_num = 0;
    std::string chain_name = "PEPTIDE1";
    for (auto& atom : rwmol.atoms()) {
        atom->setProp(MONOMER_IDX1, atom->getIdx());

        // Remove everything that doesn't have a residue number
        auto res_info =
            dynamic_cast<RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
        if (res_info) {
            to_remove.push_back(atom->getIdx());
            if (res_info->getResidueNumber() > res_num) {
                res_num = res_info->getResidueNumber();
                chain_name = res_info->getChainId();
            }
        }
    }

    rwmol.beginBatchEdit();
    for (const auto& idx : to_remove) {
        rwmol.removeAtom(rwmol.getAtomWithIdx(idx));
    }
    rwmol.commitBatchEdit();
    ++res_num;
    std::vector<std::vector<int>> frag_mapping;
    auto frags =
        RDKit::MolOps::getMolFrags(rwmol, true, nullptr, &frag_mapping);
    auto frag_idx = 0;
    for (auto frag : frags) {
        // all the atoms in this fragment (except hydrogens) belong to a new
        // residue
        for (auto at : frag->atoms()) {
            if (at->getAtomicNum() == 1) {
                continue;
            }
            auto rwmol_atom =
                rwmol.getAtomWithIdx(frag_mapping[frag_idx][at->getIdx()]);
            auto atom = mol.getAtomWithIdx(
                rwmol_atom->getProp<unsigned int>(MONOMER_IDX1));
            auto res_info = new RDKit::AtomPDBResidueInfo();
            res_info->setResidueName("X"); // this will become a SMILES monomer
            res_info->setResidueNumber(res_num);
            res_info->setChainId(chain_name);
            atom->setMonomerInfo(res_info);
        }
        ++frag_idx;
        ++res_num;
    }
}

std::pair<unsigned int, unsigned int> findFurthestResidues(
    const std::map<unsigned int, std::set<unsigned int>>& residue_graph)
{
    // Helper function to do BFS and find the furthest residue from start
    auto bfs_furthest = [&residue_graph](unsigned int start) -> unsigned int {
        std::queue<unsigned int> q;
        std::map<unsigned int, int> distance;

        q.push(start);
        distance[start] = 0;

        unsigned int furthest_res = start;
        int max_distance = 0;

        while (!q.empty()) {
            unsigned int res = q.front();
            q.pop();

            for (unsigned int neighbor : residue_graph.at(res)) {
                if (distance.find(neighbor) == distance.end()) {
                    distance[neighbor] = distance[res] + 1;
                    q.push(neighbor);

                    if (distance[neighbor] > max_distance) {
                        max_distance = distance[neighbor];
                        furthest_res = neighbor;
                    }
                }
            }
        }
        return furthest_res;
    };

    // Find two furthest residues in graph by performing a BFS twice
    // Note: this works because the graph has no cycles
    unsigned int start_node = residue_graph.begin()->first;
    auto endpoint1 = bfs_furthest(start_node);
    auto endpoint2 = bfs_furthest(endpoint1);
    return {endpoint1, endpoint2};
}

void orderResidues(RDKit::ROMol& mol)
{
    // Reorder residues by connectivity. This is needed to ensure that residue
    // attachment points are placed in the correct order.
    // Find all residues and map them to atoms
    std::map<unsigned int, std::vector<RDKit::Atom*>> residue_to_atoms;
    for (const auto& atom : mol.atoms()) {
        const auto res_info = dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
            atom->getMonomerInfo());
        if (res_info) {
            residue_to_atoms[res_info->getResidueNumber()].push_back(atom);
        }
    }

    // Create a residue connectivity graph
    std::map<unsigned int, std::set<unsigned int>> residue_connections;
    for (const auto& [res_num, atoms] : residue_to_atoms) {
        for (const auto& atom : atoms) {
            for (const auto& nbr : mol.atomNeighbors(atom)) {
                const auto nbr_res_info =
                    dynamic_cast<const RDKit::AtomPDBResidueInfo*>(
                        nbr->getMonomerInfo());
                if (nbr_res_info &&
                    static_cast<unsigned int>(
                        nbr_res_info->getResidueNumber()) != res_num) {
                    residue_connections[res_num].insert(
                        nbr_res_info->getResidueNumber());
                    residue_connections[nbr_res_info->getResidueNumber()]
                        .insert(res_num);
                }
            }
        }
    }

    // If the graph is cyclic, use input order
    // TODO: Deal with cycle and splitting into separate chains
    size_t intra_residue_bonds = 0;
    for (auto& [res_num, neighs] : residue_connections) {
        intra_residue_bonds += neighs.size();
    }
    intra_residue_bonds /= 2;
    if (intra_residue_bonds >= residue_to_atoms.size()) {
        return;
    }

    // Find the two furthest residues in the residue graph
    auto [endpoint1, endpoint2] = findFurthestResidues(residue_connections);

    // Find the lowest atom index in each endpoint residue
    unsigned int min_idx1 = std::numeric_limits<unsigned int>::max();
    unsigned int min_idx2 = std::numeric_limits<unsigned int>::max();

    for (const auto& atom : residue_to_atoms[endpoint1]) {
        min_idx1 = std::min(min_idx1, atom->getIdx());
    }

    for (const auto& atom : residue_to_atoms[endpoint2]) {
        min_idx2 = std::min(min_idx2, atom->getIdx());
    }

    // Choose the endpoint with the lowest atom index
    unsigned int start_res_num = (min_idx1 < min_idx2) ? endpoint1 : endpoint2;

    // Collect all residue numbers for array sizing
    std::set<unsigned int> residues;
    for (const auto& [res_num, _] : residue_to_atoms) {
        residues.insert(res_num);
    }
    const auto res_count = std::ranges::max(residues) + 1; // 1-based

    // Initialize data structures for BFS on residue graph
    std::vector<bool> visited_residues(res_count, false);
    std::vector<unsigned int> residue_order(res_count, 0);
    unsigned int new_res_num = 1;

    // Start BFS directly from the starting residue
    std::queue<unsigned int> residue_queue;
    residue_queue.push(start_res_num);

    visited_residues[start_res_num] = true;
    residue_order[start_res_num] = new_res_num++;

    // Perform BFS traversal of the residue graph
    while (!residue_queue.empty()) {
        unsigned int current_res = residue_queue.front();
        residue_queue.pop();

        // Process all neighboring residues
        if (residue_connections.find(current_res) !=
            residue_connections.end()) {
            for (unsigned int neighbor_res : residue_connections[current_res]) {
                if (!visited_residues[neighbor_res]) {
                    visited_residues[neighbor_res] = true;
                    residue_order[neighbor_res] = new_res_num++;
                    residue_queue.push(neighbor_res);
                }
            }
        }
    }

    // Now go through and set the new residue numbers
    for (const auto& atom : mol.atoms()) {
        const auto res_info =
            dynamic_cast<RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
        if (res_info) {
            auto old_res_num = res_info->getResidueNumber();
            res_info->setResidueNumber(residue_order[old_res_num]);
        }
    }
}

} // unnamed namespace

bool processSupGroups(RDKit::ROMol& mol)
{
    // Skip processing if there are no SUP groups
    if (!hasSupSgroups(mol)) {
        return false;
    }

    // Process the molecule with all SUP-related utilities in the right order
    removeHydrogenSupGroups(mol);
    supGroupsToPdbInfo(mol);
    addUnmappedAtomsToResidues(mol);
    orderResidues(mol);

    return true;
}

void createSupGroupsFromResidueInfo(RDKit::RWMol& mol)
{
    // Map to track (chain_id, residue_number) -> sgroup index
    std::map<std::pair<std::string, int>, unsigned int> residue_to_sgroup;
    std::vector<RDKit::SubstanceGroup> sgroups;

    // DB connection; Monomer info only stores the PDB code
    auto& monomer_db = MonomerDatabase::instance();

    // Create SUP groups for each unique residue
    for (const auto& at : mol.atoms()) {
        auto* monomer_info =
            dynamic_cast<RDKit::AtomPDBResidueInfo*>(at->getMonomerInfo());
        if (!monomer_info) {
            continue; // Skip atoms without PDB residue info
        }

        std::string chain_id = monomer_info->getChainId();
        // For now, only PEPTIDE monomers are written to SUP groups
        if (chain_id.find("PEPTIDE") != 0) {
            continue;
        }

        int residue_number = monomer_info->getResidueNumber();
        std::string residue_name = monomer_info->getResidueName();

        auto helm_info = monomer_db.getHelmInfo(residue_name);
        std::string monomer_label = "X";
        if (helm_info) {
            monomer_label = std::get<0>(*helm_info);
        }

        // Chains are not identified by anything in SUP groups, we just need
        // to make sure monomers from different chains get different SUP groups
        auto residue_key = std::make_pair(chain_id, residue_number);

        // Check if we've already created a SUP group for this residue
        if (residue_to_sgroup.find(residue_key) == residue_to_sgroup.end()) {
            // Create a new SUP group for this residue
            RDKit::SubstanceGroup sgroup(&mol, "SUP");
            sgroup.setProp(SGROUP_PROP_CLASS, "AA");
            sgroup.setProp(SGROUP_PROP_LABEL, monomer_label);
            sgroup.setProp<unsigned int>(
                SGROUP_PROP_RESNUM,
                static_cast<unsigned int>(sgroups.size() + 1));

            unsigned int sgroup_idx = sgroups.size();
            sgroups.push_back(sgroup);
            residue_to_sgroup[residue_key] = sgroup_idx;
        }

        // Add this atom to its residue's SUP group
        unsigned int sgroup_idx = residue_to_sgroup[residue_key];
        sgroups[sgroup_idx].addAtomWithIdx(at->getIdx());
    }

    // Add all the SUP groups to the molecule
    auto& mol_sgroups = RDKit::getSubstanceGroups(mol);
    mol_sgroups = std::move(sgroups);
}

} // namespace rdkit_extensions
} // namespace schrodinger