#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

using AnalogMap =
    std::unordered_map<std::string, std::vector<rdkit_extensions::MonomerInfo>>;

/**
 * One unit of work for `MolModel::mutateMonomers`: a set of atoms that
 * should all be mutated to the same target HELM symbol.
 */
struct MonomerMutation {
    std::unordered_set<const RDKit::Atom*> atoms;
    std::string helm_symbol;
};

/**
 * True iff `atoms` is non-empty and every atom has the given monomer type.
 * Empty returns false (gates menu actions: nothing-to-act-on must be false).
 */
SKETCHER_API bool
all_monomers_have_type(const std::unordered_set<const RDKit::Atom*>& atoms,
                       MonomerType type);

/**
 * True iff every atom in `atoms` is currently a D-form residue.
 */
SKETCHER_API bool
all_residues_are_d_form(const std::unordered_set<const RDKit::Atom*>& atoms);

/**
 * `dFoo` is the D-form of `Foo` iff both exist as PEPTIDEs and share the
 * same NATURAL_ANALOG. Prefix alone is unsafe — a custom DB could name
 * an entry `dXyz` with no semantic link to `Xyz`.
 */
SKETCHER_API bool is_d_form(std::string_view sym);

/**
 * Sugar parallel of `is_d_form`: `dX` is the deoxy form of `X` iff both
 * exist in the DB and share the same NATURAL_ANALOG.
 */
SKETCHER_API bool is_deoxy_sugar(std::string_view sym);

/**
 * The L↔D toggle target for `sym`, or nullopt when no DB-corroborated
 * counterpart exists (would silently change chemistry otherwise).
 */
SKETCHER_API std::optional<std::string>
toggle_d_form_symbol(std::string_view sym);

/**
 * The R↔dR toggle target for `sym`, or nullopt when the candidate is
 * missing or has a different NATURAL_ANALOG.
 */
SKETCHER_API std::optional<std::string>
toggle_sugar_symbol(std::string_view sym);

/**
 * Build the per-target-symbol mutation batch for the D-form toggle.
 * Atoms already at the target form are skipped; atoms whose toggle has
 * no DB counterpart are dropped. One MonomerMutation per distinct target
 * symbol.
 */
SKETCHER_API std::vector<MonomerMutation>
build_d_form_batch(const std::unordered_set<const RDKit::Atom*>& atoms,
                   bool target_d_form);

/**
 * Sugar-toggle batch with direction fixed by the caller. Skips atoms
 * already at `target_sym` and atoms whose toggle leads elsewhere.
 */
SKETCHER_API std::vector<MonomerMutation>
build_sugar_toggle_batch(const std::unordered_set<const RDKit::Atom*>& atoms,
                         const std::string& target_sym);

/**
 * Filter to the atoms whose current symbol differs from `target` — avoids
 * pushing a no-op mutation (and a useless undo entry).
 */
SKETCHER_API std::unordered_set<const RDKit::Atom*>
atoms_needing_mutation(const std::unordered_set<const RDKit::Atom*>& atoms,
                       std::string_view target);

/**
 * NA analogs keyed by natural-analog symbol, merged across RNA and DNA
 * chain types and deduped by analog symbol — custom DBs may register
 * variants under either polymer type.
 */
SKETCHER_API AnalogMap get_merged_na_analogs();

/**
 * Build the analog list for an NA button: optionally prepend a parent
 * standard (e.g. R for the dR button) sourced from `standard_names`,
 * then append the natural-analog list keyed by the override target (or
 * `symbol` itself if no override applies), filtering out entries whose
 * symbol matches `symbol` (the button's own standard).
 */
SKETCHER_API std::vector<rdkit_extensions::MonomerInfo>
get_analogs_for_na_button(
    const std::string& symbol, const AnalogMap& analogs_by_na,
    const std::unordered_map<std::string, std::string>& standard_names);

} // namespace sketcher
} // namespace schrodinger
