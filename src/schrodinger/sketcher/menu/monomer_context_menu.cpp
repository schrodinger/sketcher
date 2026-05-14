#include "schrodinger/sketcher/menu/monomer_context_menu.h"

#include <algorithm>
#include <string>

#include <QAction>
#include <QMenu>
#include <rdkit/GraphMol/Atom.h>

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/monomer_analog.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

// True when every selected monomer is part of a nucleic acid chain (no
// peptides, no CHEM) AND at least one is an NA base — the gate for the
// "Add Complementary Sequence" action visibility.
static bool na_selection_with_at_least_one_base(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    // Any non-NA monomer type (peptide, CHEM, or an unknown future variant)
    // disqualifies the whole selection.
    auto is_na_monomer = [](const RDKit::Atom* a) {
        switch (get_monomer_type(a)) {
            case MonomerType::NA_BASE:
            case MonomerType::NA_SUGAR:
            case MonomerType::NA_PHOSPHATE:
                return true;
            default:
                return false;
        }
    };
    auto is_na_base = [](const RDKit::Atom* a) {
        return get_monomer_type(a) == MonomerType::NA_BASE;
    };
    // any_of is false for an empty selection, so no explicit empty guard.
    return std::ranges::all_of(atoms, is_na_monomer) &&
           std::ranges::any_of(atoms, is_na_base);
}

// Filter the menu's selection down to NA bases whose symbol has a
// Watson-Crick complement. Used to drive both the "enabled" state of
// the action and the payload of the emitted signal.
static std::unordered_set<const RDKit::Atom*>
complementable_bases(const std::unordered_set<const RDKit::Atom*>& atoms)
{
    std::unordered_set<const RDKit::Atom*> out;
    for (const auto* a : atoms) {
        if (get_monomer_type(a) != MonomerType::NA_BASE) {
            continue;
        }
        if (na_base_has_complement(get_monomer_res_name(a))) {
            out.insert(a);
        }
    }
    return out;
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
    createAddComplementaryStrandAction();
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
        emit mutateMonomerRequested({{std::move(target_atoms), std::move(sym)}},
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
        emit mutateMonomerRequested(
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
    auto analogs_by_base = get_merged_na_analogs();

    auto emit_mutation = [this](std::string sym) {
        auto target_atoms = atoms_needing_mutation(m_atoms, sym);
        if (target_atoms.empty()) {
            return;
        }
        emit mutateMonomerRequested({{std::move(target_atoms), std::move(sym)}},
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
        const auto description = QString("Change %1 to %2")
                                     .arg(QString::fromStdString(from))
                                     .arg(QString::fromStdString(*target));
        emit mutateMonomerRequested(std::move(batch), description);
    });
}

void MonomerContextMenu::createAddComplementaryStrandAction()
{
    m_add_complement_action =
        addAction("Add Complementary Sequence", this, [this]() {
            // Reuse the set computed in updateActions() rather than walking
            // the selection a third time.
            if (m_complement_bases.empty()) {
                return;
            }
            emit addComplementaryStrandRequested(m_complement_bases);
        });
}

void MonomerContextMenu::createDeleteAction()
{
    addAction("Delete", this, [this]() { emit deleteRequested(m_atoms); });
}

void MonomerContextMenu::updateActions()
{
    const bool all_peptide =
        all_monomers_have_type(m_atoms, MonomerType::PEPTIDE);
    const bool all_na_base =
        all_monomers_have_type(m_atoms, MonomerType::NA_BASE);
    const bool all_na_sugar =
        all_monomers_have_type(m_atoms, MonomerType::NA_SUGAR);

    m_mutate_residue_menu->menuAction()->setVisible(all_peptide);
    m_set_d_form_action->setVisible(all_peptide);
    m_protonate_action->setVisible(all_peptide);
    m_mutate_base_menu->menuAction()->setVisible(all_na_base);
    m_sugar_toggle_action->setVisible(all_na_sugar);

    m_complement_bases = complementable_bases(m_atoms);
    const bool na_with_base = na_selection_with_at_least_one_base(m_atoms);
    m_add_complement_action->setVisible(na_with_base);
    // Enabled iff at least one selected base has a Watson-Crick complement
    // symbol; DB-level validation happens model-side.
    m_add_complement_action->setEnabled(na_with_base &&
                                        !m_complement_bases.empty());

    if (all_peptide) {
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

    if (all_na_sugar) {
        // Direction comes from the right-clicked atom so a mixed {R, dR}
        // selection acts on whatever the user clicked. Disabled when the
        // primary atom is unavailable — no unambiguous direction.
        std::optional<std::string> target;
        std::string from;
        if (m_primary_atom != nullptr &&
            get_monomer_type(m_primary_atom) == MonomerType::NA_SUGAR) {
            from = get_monomer_res_name(m_primary_atom);
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
