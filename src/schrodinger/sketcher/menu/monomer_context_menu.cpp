#include "schrodinger/sketcher/menu/monomer_context_menu.h"

#include <algorithm>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>

#include <QAction>
#include <QMenu>
#include <rdkit/GraphMol/Atom.h>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

static bool all_peptide(const std::unordered_set<const RDKit::Atom*>& atoms)
{
    if (atoms.empty()) {
        return false;
    }
    return std::ranges::all_of(atoms, [](const auto* a) {
        return get_monomer_type(a) == MonomerType::PEPTIDE;
    });
}

static bool
all_of_monomer_type(const std::unordered_set<const RDKit::Atom*>& atoms,
                    MonomerType type)
{
    if (atoms.empty()) {
        return false;
    }
    return std::ranges::all_of(
        atoms, [type](const auto* a) { return get_monomer_type(a) == type; });
}

static rdkit_extensions::opt_string_t
get_peptide_natural_analog(std::string_view sym)
{
    return rdkit_extensions::MonomerDatabase::instance().getNaturalAnalog(
        std::string(sym), rdkit_extensions::ChainType::PEPTIDE);
}

// `dFoo` is the D-form of `Foo` iff both exist as PEPTIDEs and share
// the same NATURAL_ANALOG. Prefix alone is unsafe — a custom DB could
// name an entry `dXyz` with no semantic link to `Xyz`.
static bool is_d_form(std::string_view sym)
{
    if (sym.size() < 2 || sym.front() != 'd') {
        return false;
    }
    auto sym_analog = get_peptide_natural_analog(sym);
    auto stripped_analog = get_peptide_natural_analog(sym.substr(1));
    return sym_analog && stripped_analog && !sym_analog->empty() &&
           *sym_analog == *stripped_analog;
}

static std::optional<std::string> toggle_d_form_symbol(std::string_view sym)
{
    if (is_d_form(sym)) {
        // Stripping the prefix is safe: is_d_form only succeeds when
        // the L-form counterpart exists with a matching analog.
        return std::string(sym.substr(1));
    }
    // L → D: candidate "d" + sym is a true D-form when both sides
    // exist in the DB AND share the same NATURAL_ANALOG.
    auto candidate = std::string("d") + std::string(sym);
    auto candidate_analog = get_peptide_natural_analog(candidate);
    auto sym_analog = get_peptide_natural_analog(sym);
    if (candidate_analog && sym_analog && !candidate_analog->empty() &&
        *candidate_analog == *sym_analog) {
        return candidate;
    }
    return std::nullopt;
}

// True only when every atom in `atoms` is already a D-form residue.
// Used to drive the Set D-/L-Form action label and target: when this
// returns true, the action reads "Set L-Form" and clicking it flips
// the residues back to L; otherwise it reads "Set D-Form".
static bool
all_residues_are_d_form(const std::unordered_set<const RDKit::Atom*>& atoms)
{
    return std::ranges::all_of(atoms, [](const auto* a) {
        return is_d_form(get_monomer_res_name(a));
    });
}

// Build the per-target-symbol mutation batch for the D-form toggle.
// Atoms already at the target form are skipped; atoms whose toggle has
// no DB counterpart are dropped. Returns one MonomerMutation per
// distinct target symbol so the receiver can coalesce them under one
// undo entry.
static std::vector<MonomerMutation>
build_d_form_batch(const std::unordered_set<const RDKit::Atom*>& atoms,
                   bool target_d_form)
{
    std::unordered_map<std::string, std::unordered_set<const RDKit::Atom*>>
        by_target;
    for (const auto* a : atoms) {
        const auto current = get_monomer_res_name(a);
        if (is_d_form(current) == target_d_form) {
            continue;
        }
        if (auto next = toggle_d_form_symbol(current)) {
            by_target[*next].insert(a);
        }
    }
    std::vector<MonomerMutation> batch;
    batch.reserve(by_target.size());
    for (auto& [sym, group] : by_target) {
        batch.push_back({std::move(group), sym});
    }
    return batch;
}

// Filter to the atoms whose current symbol differs from `target` —
// avoids pushing a no-op mutation (and a useless undo entry) when the
// user picks a target that the residue is already at.
static std::unordered_set<const RDKit::Atom*>
atoms_needing_mutation(const std::unordered_set<const RDKit::Atom*>& atoms,
                       std::string_view target)
{
    std::unordered_set<const RDKit::Atom*> out;
    for (const auto* a : atoms) {
        if (get_monomer_res_name(a) != target) {
            out.insert(a);
        }
    }
    return out;
}

// Default DB tags sugars as RNA, but custom DBs may register them under
// DNA — try RNA first, then fall back to DNA so either works.
static rdkit_extensions::opt_string_t
get_na_natural_analog(std::string_view sym)
{
    auto& db = rdkit_extensions::MonomerDatabase::instance();
    if (auto rna = db.getNaturalAnalog(std::string(sym),
                                       rdkit_extensions::ChainType::RNA);
        rna && !rna->empty()) {
        return rna;
    }
    return db.getNaturalAnalog(std::string(sym),
                               rdkit_extensions::ChainType::DNA);
}

// Sugar parallel of `is_d_form`: `dX` is the deoxy form of `X` iff both
// exist in the DB and share the same NATURAL_ANALOG. Prefix alone is
// unsafe — a custom `dXyz` might be unrelated to `Xyz`.
static bool is_deoxy_sugar(std::string_view sym)
{
    if (sym.size() < 2 || sym.front() != 'd') {
        return false;
    }
    auto sym_analog = get_na_natural_analog(sym);
    auto stripped_analog = get_na_natural_analog(sym.substr(1));
    return sym_analog && stripped_analog && !sym_analog->empty() &&
           *sym_analog == *stripped_analog;
}

// Returns the symbol on the other side of the R↔dR toggle for `sym`,
// or nullopt when no DB-corroborated counterpart exists.
static std::optional<std::string> toggle_sugar_symbol(std::string_view sym)
{
    if (is_deoxy_sugar(sym)) {
        return std::string(sym.substr(1));
    }
    auto candidate = std::string("d") + std::string(sym);
    auto candidate_analog = get_na_natural_analog(candidate);
    auto sym_analog = get_na_natural_analog(sym);
    if (candidate_analog && sym_analog && !candidate_analog->empty() &&
        *candidate_analog == *sym_analog) {
        return candidate;
    }
    return std::nullopt;
}

// Sugar-toggle batch with direction fixed by the caller. Skips atoms
// already at `target_sym` and atoms whose toggle leads elsewhere
// (defensive against custom DBs that mix unrelated sugar symbols).
static std::vector<MonomerMutation>
build_sugar_toggle_batch(const std::unordered_set<const RDKit::Atom*>& atoms,
                         const std::string& target_sym)
{
    std::unordered_set<const RDKit::Atom*> group;
    for (const auto* a : atoms) {
        if (get_monomer_res_name(a) == target_sym) {
            continue;
        }
        if (auto next = toggle_sugar_symbol(get_monomer_res_name(a));
            next && *next == target_sym) {
            group.insert(a);
        }
    }
    if (group.empty()) {
        return {};
    }
    return {{std::move(group), target_sym}};
}

// Base analogs merged across RNA and DNA chains, deduped by symbol —
// custom DBs may register variants under either polymer type. Same
// pattern as `monomer_tool_widget.cpp::get_merged_na_analogs`.
static std::unordered_map<std::string,
                          std::vector<rdkit_extensions::MonomerInfo>>
get_merged_base_analogs()
{
    auto rna =
        rdkit_extensions::MonomerDatabase::instance()
            .getMonomersByNaturalAnalog(rdkit_extensions::ChainType::RNA);
    auto dna =
        rdkit_extensions::MonomerDatabase::instance()
            .getMonomersByNaturalAnalog(rdkit_extensions::ChainType::DNA);
    for (auto& [analog, infos] : dna) {
        auto& dest = rna[analog];
        std::unordered_set<std::string> seen;
        for (const auto& info : dest) {
            if (info.symbol) {
                seen.insert(*info.symbol);
            }
        }
        for (auto& info : infos) {
            if (info.symbol && !seen.contains(*info.symbol)) {
                seen.insert(*info.symbol);
                dest.push_back(std::move(info));
            }
        }
    }
    return rna;
}

MonomerContextMenu::MonomerContextMenu(QWidget* parent) :
    AbstractContextMenu(parent)
{
    setTitle("Monomer");
    createMutateResidueSubMenu();
    createSetDFormAction();
    createProtonateAction();
    createMutateBaseSubMenu();
    createSugarToggleAction();
    createDeleteAction();
}

void MonomerContextMenu::createMutateResidueSubMenu()
{
    // Palette of natural AAs each with a sub-submenu of analogs queried
    // from the monomer DB. Captured once at construction; mid-session
    // custom-DB updates won't reflect until the menu is rebuilt.
    m_mutate_residue_menu = addMenu("Mutate Residue");
    auto analogs_by_aa =
        rdkit_extensions::MonomerDatabase::instance()
            .getMonomersByNaturalAnalog(rdkit_extensions::ChainType::PEPTIDE);

    auto emit_mutation = [this](std::string sym) {
        auto target_atoms = atoms_needing_mutation(m_atoms, sym);
        if (target_atoms.empty()) {
            return;
        }
        emit mutateResidueRequested({{std::move(target_atoms), std::move(sym)}},
                                    {});
    };

    for (auto aa_tool : AMINO_ACID_TOOL_DISPLAY_ORDER) {
        const auto& symbol = AMINO_ACID_TOOL_TO_RES_NAME.at(aa_tool);
        const auto& name = AMINO_ACID_TOOL_TO_FULL_NAME.at(aa_tool);
        const auto title = QString::fromStdString(name + " (" + symbol + ")");
        auto* sub = m_mutate_residue_menu->addMenu(title);

        sub->addAction(QString::fromStdString(symbol), this,
                       [emit_mutation, sym = symbol]() { emit_mutation(sym); });

        const auto it = analogs_by_aa.find(symbol);
        if (it == analogs_by_aa.end() || it->second.empty()) {
            continue;
        }
        auto analogs = it->second;
        std::ranges::sort(analogs, {}, [](const auto& info) {
            return info.symbol.value_or("");
        });
        for (const auto& info : analogs) {
            if (!info.symbol) {
                continue;
            }
            const auto sym = *info.symbol;
            QString label = QString::fromStdString(sym);
            if (info.name && !info.name->empty()) {
                label +=
                    QString(" — %1").arg(QString::fromStdString(*info.name));
            }
            sub->addAction(label, this,
                           [emit_mutation, sym]() { emit_mutation(sym); });
        }
    }
}

void MonomerContextMenu::createSetDFormAction()
{
    m_set_d_form_action = addAction("Set D-Form", this, [this]() {
        const bool target_d_form = !all_residues_are_d_form(m_atoms);
        auto batch = build_d_form_batch(m_atoms, target_d_form);
        if (batch.empty()) {
            return;
        }
        emit mutateResidueRequested(
            std::move(batch), target_d_form ? "Set D-Form" : "Set L-Form");
    });
}

void MonomerContextMenu::createProtonateAction()
{
    // Placeholder for a future SKETCH-2648 phase — the monomer DB has
    // no protonated/deprotonated symbol pairs today, so there's nothing
    // to wire. Shown disabled to reserve the slot per the spec mockup.
    m_protonate_action = addAction("Protonate");
    m_protonate_action->setEnabled(false);
}

void MonomerContextMenu::createMutateBaseSubMenu()
{
    // Same shape as Mutate Residue. Analogs captured once at construction;
    // mid-session custom-DB updates won't reflect until the menu rebuilds.
    m_mutate_base_menu = addMenu("Mutate Base");
    auto analogs_by_base = get_merged_base_analogs();

    auto emit_mutation = [this](std::string sym) {
        auto target_atoms = atoms_needing_mutation(m_atoms, sym);
        if (target_atoms.empty()) {
            return;
        }
        emit mutateResidueRequested({{std::move(target_atoms), std::move(sym)}},
                                    {});
    };

    for (auto base_tool : NA_BASE_DISPLAY_ORDER) {
        const auto& symbol = NUCLEIC_ACID_TOOL_TO_RES_NAME.at(base_tool);
        const auto& name = NA_BASE_TO_FULL_NAME.at(base_tool);
        const auto title = QString::fromStdString(name + " (" + symbol + ")");
        auto* sub = m_mutate_base_menu->addMenu(title);

        sub->addAction(QString::fromStdString(symbol), this,
                       [emit_mutation, sym = symbol]() { emit_mutation(sym); });

        const auto it = analogs_by_base.find(symbol);
        if (it == analogs_by_base.end() || it->second.empty()) {
            continue;
        }
        auto analogs = it->second;
        std::ranges::sort(analogs, {}, [](const auto& info) {
            return info.symbol.value_or("");
        });
        for (const auto& info : analogs) {
            if (!info.symbol) {
                continue;
            }
            const auto sym = *info.symbol;
            QString label = QString::fromStdString(sym);
            if (info.name && !info.name->empty()) {
                label +=
                    QString(" — %1").arg(QString::fromStdString(*info.name));
            }
            sub->addAction(label, this,
                           [emit_mutation, sym]() { emit_mutation(sym); });
        }
    }
}

void MonomerContextMenu::createSugarToggleAction()
{
    m_sugar_toggle_action = addAction("Change R to dR", this, [this]() {
        if (m_primary_atom == nullptr) {
            return;
        }
        const auto target =
            toggle_sugar_symbol(get_monomer_res_name(m_primary_atom));
        if (!target) {
            return;
        }
        auto batch = build_sugar_toggle_batch(m_atoms, *target);
        if (batch.empty()) {
            return;
        }
        const auto from = get_monomer_res_name(m_primary_atom);
        const auto description =
            QString("Change %1 to %2")
                .arg(QString::fromStdString(std::string(from)))
                .arg(QString::fromStdString(*target));
        emit mutateResidueRequested(std::move(batch), description);
    });
}

void MonomerContextMenu::createDeleteAction()
{
    addAction("Delete", this, [this]() { emit deleteRequested(m_atoms); });
}

void MonomerContextMenu::updateActions()
{
    const bool peptide = all_peptide(m_atoms);
    const bool na_base = all_of_monomer_type(m_atoms, MonomerType::NA_BASE);
    const bool na_sugar = all_of_monomer_type(m_atoms, MonomerType::NA_SUGAR);

    m_mutate_residue_menu->menuAction()->setVisible(peptide);
    m_set_d_form_action->setVisible(peptide);
    m_protonate_action->setVisible(peptide);
    m_mutate_base_menu->menuAction()->setVisible(na_base);
    m_sugar_toggle_action->setVisible(na_sugar);

    if (peptide) {
        bool any_d_form_toggleable = false;
        for (const auto* a : m_atoms) {
            if (toggle_d_form_symbol(get_monomer_res_name(a)).has_value()) {
                any_d_form_toggleable = true;
                break;
            }
        }
        const bool target_d_form = !all_residues_are_d_form(m_atoms);
        m_set_d_form_action->setText(target_d_form ? "Set D-Form"
                                                   : "Set L-Form");
        m_set_d_form_action->setEnabled(any_d_form_toggleable);
    }

    if (na_sugar) {
        // Direction comes from the right-clicked atom so a mixed {R, dR}
        // selection acts on whatever the user clicked. Disabled when the
        // primary atom is unavailable — no unambiguous direction.
        std::optional<std::string> target;
        std::string from;
        if (m_primary_atom != nullptr &&
            get_monomer_type(m_primary_atom) == MonomerType::NA_SUGAR) {
            from = std::string(get_monomer_res_name(m_primary_atom));
            target = toggle_sugar_symbol(from);
        }
        if (target) {
            m_sugar_toggle_action->setText(
                QString("Change %1 to %2")
                    .arg(QString::fromStdString(from))
                    .arg(QString::fromStdString(*target)));
            m_sugar_toggle_action->setEnabled(
                !build_sugar_toggle_batch(m_atoms, *target).empty());
        } else {
            m_sugar_toggle_action->setText("Change R to dR");
            m_sugar_toggle_action->setEnabled(false);
        }
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/monomer_context_menu.moc"
