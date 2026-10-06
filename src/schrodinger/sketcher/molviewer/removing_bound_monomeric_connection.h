#pragma once

#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include <QString>

#include "schrodinger/sketcher/definitions.h"

class QWidget;

namespace RDKit
{
class Atom;
class ROMol;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class MessageBoxDialog;

/**
 * An attachment point that is involved in a bound connection.  If a mutation
 * removes this attachment point, then the user should be warned that the
 * associated connections will be automatically removed.
 */
struct RequiredConnection {
    int attachment_point;
    unsigned int bound_monomer_index;
    bool is_secondary_connection;
};

/**
 * Atom indices and the HELM symbol to mutate them to.  Note that atom indices
 * remain valid across snapshot-based monomer mutations.
 */
struct IndexedMonomerMutation {
    std::vector<unsigned int> atom_indices;
    std::string helm_symbol;
};

/**
 * Connections that a mutation would remove due to removing the bound
 * attachment point on the mutated monomer.
 */
struct MutationConnectionChanges {
    std::vector<int> missing_attachment_points;
    std::unordered_set<unsigned int> bond_indices;
    std::unordered_set<unsigned int> secondary_connection_indices;
    unsigned int affected_monomers = 0;
};

/**
 * Return a monomer's bound numbered attachment points and their connections.
 */
SKETCHER_API std::pair<std::vector<int>, std::vector<RequiredConnection>>
get_required_attachment_points(const RDKit::Atom* monomer);

/**
 * Find connections that would lose their attachment points after mutation.
 *
 * @param mol The current monomeric molecule, before any mutation or bond
 * removal. Every atom index in `mutations` must refer to an atom in this mol.
 * @param mutations Groups of atom indices and replacement HELM symbols. Each
 * symbol is looked up with its atom's chain type to determine which numbered
 * attachment points the replacement provides.
 * @return Missing bound attachment point numbers (sorted and unique), indices
 * of primary and secondary connections to remove, and the number of affected
 * monomers. Bond indices are captured before mutation so the caller can
 * resolve them to live bonds immediately before removal.
 */
SKETCHER_API MutationConnectionChanges get_monomer_mutation_connection_changes(
    const RDKit::ROMol& mol,
    const std::vector<IndexedMonomerMutation>& mutations);

/**
 * Format numbered R groups for monomer attachment-point messages.
 */
SKETCHER_API QString
format_r_group_list(const std::vector<int>& r_group_numbers);

/**
 * Show the standard warning before removing connections at bound R groups.
 */
SKETCHER_API MessageBoxDialog*
show_bound_connection_warning(const std::vector<int>& missing_attachment_points,
                              QWidget* parent, bool multiple_monomers = false);

} // namespace sketcher
} // namespace schrodinger
