#include "schrodinger/sketcher/molviewer/removing_bound_monomeric_connection.h"

#include <algorithm>

#include <QStringList>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/dialog/message_box_dialog.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

std::pair<std::vector<int>, std::vector<RequiredConnection>>
get_required_attachment_points(const RDKit::Atom* monomer)
{
    std::vector<int> required_attachment_points;
    std::vector<RequiredConnection> required_connections;
    const auto attachment_points = get_attachment_points_for_monomer(monomer);
    for (const auto& bound_attachment_point : attachment_points.first) {
        if (bound_attachment_point.num <= 0) {
            continue;
        }
        required_attachment_points.push_back(bound_attachment_point.num);
        required_connections.push_back(
            {bound_attachment_point.num,
             bound_attachment_point.bound_monomer->getIdx(),
             bound_attachment_point.is_secondary_connection});
    }
    return std::make_pair(required_attachment_points, required_connections);
}

MutationConnectionChanges get_monomer_mutation_connection_changes(
    const RDKit::ROMol& mol,
    const std::vector<IndexedMonomerMutation>& mutations)
{
    MutationConnectionChanges changes;
    for (const auto& [atom_indices, symbol] : mutations) {
        for (const auto atom_index : atom_indices) {
            const auto* atom = mol.getAtomWithIdx(atom_index);
            const auto [required_points, connections] =
                get_required_attachment_points(atom);
            if (required_points.empty()) {
                continue;
            }
            const auto available_points = get_attachment_points_for_res(
                symbol, rdkit_extensions::getChainType(*atom));
            std::unordered_set<int> available_numbers;
            for (const auto& [number, element] : available_points) {
                available_numbers.insert(number);
            }
            bool missing_point_on_monomer = false;
            for (const auto& connection : connections) {
                if (available_numbers.contains(connection.attachment_point)) {
                    continue;
                }
                changes.missing_attachment_points.push_back(
                    connection.attachment_point);
                missing_point_on_monomer = true;
                const auto* bond = mol.getBondBetweenAtoms(
                    atom_index, connection.bound_monomer_index);
                auto& indices = connection.is_secondary_connection
                                    ? changes.secondary_connection_indices
                                    : changes.bond_indices;
                indices.insert(bond->getIdx());
            }
            if (missing_point_on_monomer) {
                ++changes.affected_monomers;
            }
        }
    }
    std::ranges::sort(changes.missing_attachment_points);
    changes.missing_attachment_points.erase(
        std::unique(changes.missing_attachment_points.begin(),
                    changes.missing_attachment_points.end()),
        changes.missing_attachment_points.end());
    return changes;
}

QString format_r_group_list(const std::vector<int>& r_group_numbers)
{
    if (r_group_numbers.empty()) {
        return {};
    }

    QStringList r_groups;
    for (const auto r_group_num : r_group_numbers) {
        r_groups.append("R" + QString::number(r_group_num));
    }

    if (r_groups.size() == 1) {
        return r_groups.front();
    }

    const auto last_r_group = r_groups.takeLast();
    const auto separator = r_groups.size() == 1 ? " " : ", ";
    return r_groups.join(", ") + separator + "and " + last_r_group;
}

MessageBoxDialog*
show_bound_connection_warning(const std::vector<int>& missing_attachment_points,
                              QWidget* parent, bool multiple_monomers)
{
    const bool plural = missing_attachment_points.size() != 1;
    const auto attachment_points =
        format_r_group_list(missing_attachment_points);
    const auto warning_text =
        multiple_monomers
            ? "The selected replacement monomers are missing bound " +
                  attachment_points +
                  " attachment points. Continuing will remove the affected "
                  "connections."
            : attachment_points + (plural ? " have" : " has") +
                  " been removed from this monomer but " +
                  (plural ? "are" : "is") +
                  " currently bound. Continuing will remove " +
                  (plural ? "these connections." : "this connection.");
    return show_warning_dialog("Remove Bound Connections?", warning_text,
                               parent);
}

} // namespace sketcher
} // namespace schrodinger
