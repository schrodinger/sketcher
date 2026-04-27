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

namespace
{

bool all_peptide(const std::unordered_set<const RDKit::Atom*>& atoms)
{
    if (atoms.empty()) {
        return false;
    }
    return std::ranges::all_of(atoms, [](const auto* a) {
        return get_monomer_type(a) == MonomerType::PEPTIDE;
    });
}

rdkit_extensions::opt_string_t peptide_natural_analog(std::string_view sym)
{
    return rdkit_extensions::MonomerDatabase::instance().getNaturalAnalog(
        std::string(sym), rdkit_extensions::ChainType::PEPTIDE);
}

// `dFoo` is the D-form of `Foo` iff both exist as PEPTIDEs and share
// the same NATURAL_ANALOG. Prefix alone is unsafe — a custom DB could
// name an entry `dXyz` with no semantic link to `Xyz`.
bool is_d_form(std::string_view sym)
{
    if (sym.size() < 2 || sym.front() != 'd') {
        return false;
    }
    auto sym_analog = peptide_natural_analog(sym);
    auto stripped_analog = peptide_natural_analog(sym.substr(1));
    return sym_analog && stripped_analog && !sym_analog->empty() &&
           *sym_analog == *stripped_analog;
}

std::optional<std::string> toggle_d_form_symbol(std::string_view sym)
{
    if (is_d_form(sym)) {
        // Stripping the prefix is safe: is_d_form only succeeds when
        // the L-form counterpart exists with a matching analog.
        return std::string(sym.substr(1));
    }
    // L → D: candidate "d" + sym is a true D-form when both sides
    // exist in the DB AND share the same NATURAL_ANALOG.
    auto candidate = std::string("d") + std::string(sym);
    auto candidate_analog = peptide_natural_analog(candidate);
    auto sym_analog = peptide_natural_analog(sym);
    if (candidate_analog && sym_analog && !candidate_analog->empty() &&
        *candidate_analog == *sym_analog) {
        return candidate;
    }
    return std::nullopt;
}

// Target D-Form unless every atom is already D-form.
bool should_target_d_form(const std::unordered_set<const RDKit::Atom*>& atoms)
{
    return !std::ranges::all_of(atoms, [](const auto* a) {
        return is_d_form(get_monomer_res_name(a));
    });
}

// Build the per-target-symbol mutation batch for the D-form toggle.
// Atoms already at the target form are skipped; atoms whose toggle has
// no DB counterpart are dropped. Returns one MonomerMutation per
// distinct target symbol so the receiver can coalesce them under one
// undo entry.
std::vector<MonomerMutation>
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
std::unordered_set<const RDKit::Atom*>
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

} // namespace

MonomerContextMenu::MonomerContextMenu(QWidget* parent) :
    AbstractContextMenu(parent)
{
    setTitle("Monomer");
    createMutateResidueSubMenu();
    createSetDFormAction();
    createProtonateAction();
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

        const auto natural_sym = symbol;
        sub->addAction(
            QString::fromStdString(natural_sym), this,
            [emit_mutation, natural_sym]() { emit_mutation(natural_sym); });

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
        const bool target_d_form = should_target_d_form(m_atoms);
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

void MonomerContextMenu::createDeleteAction()
{
    addAction("Delete", this, [this]() { emit deleteRequested(m_atoms); });
}

void MonomerContextMenu::updateActions()
{
    const bool peptide = all_peptide(m_atoms);

    m_mutate_residue_menu->menuAction()->setVisible(peptide);
    m_set_d_form_action->setVisible(peptide);
    m_protonate_action->setVisible(peptide);

    if (!peptide) {
        return;
    }

    bool any_d_form_toggleable = false;
    for (const auto* a : m_atoms) {
        if (toggle_d_form_symbol(get_monomer_res_name(a)).has_value()) {
            any_d_form_toggleable = true;
            break;
        }
    }
    const bool target_d_form = should_target_d_form(m_atoms);
    m_set_d_form_action->setText(target_d_form ? "Set D-Form" : "Set L-Form");
    m_set_d_form_action->setEnabled(any_d_form_toggleable);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/monomer_context_menu.moc"
