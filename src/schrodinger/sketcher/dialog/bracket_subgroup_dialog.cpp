#include "schrodinger/sketcher/dialog/bracket_subgroup_dialog.h"

#include <iostream>
#include <unordered_map>

#include <QString>

#include "schrodinger/sketcher/rdkit/sgroup.h"
#include "schrodinger/sketcher/numeric_label_validator.h"
#include "schrodinger/sketcher/rdkit/s_group_constants.h"
#include "schrodinger/sketcher/ui/ui_bracket_subgroup_dialog.h"

Q_DECLARE_METATYPE(schrodinger::sketcher::RepeatPattern);
Q_DECLARE_METATYPE(schrodinger::sketcher::SubgroupType);

namespace
{

using schrodinger::sketcher::RepeatPattern;
using schrodinger::sketcher::SubgroupType;

const std::vector<std::pair<SubgroupType, QString>> subgrouptype_titles = {
    {SubgroupType::SRU_POLYMER, "SRU polymer"},
    {SubgroupType::COPOLYMER, "Copolymer"},
};

const std::vector<std::pair<RepeatPattern, QString>> repeatpattern_titles = {
    {RepeatPattern::HEAD_TO_TAIL, "Head-to-tail"},
    {RepeatPattern::HEAD_TO_HEAD, "Head-to-head"},
    {RepeatPattern::EITHER_UNKNOWN, "Either / Unknown"},
};

const std::unordered_map<SubgroupType, QString> subgrouptype_label_map = {
    {SubgroupType::OTHER, ""},
    {SubgroupType::SRU_POLYMER, ""},
    {SubgroupType::COPOLYMER, "co"},
};

} // anonymous namespace

namespace schrodinger
{
namespace sketcher
{

AbstractBracketSubgroupDialog::AbstractBracketSubgroupDialog(QWidget* parent) :
    ModalDialog(parent)
{
    ui.reset(new Ui::BracketSubgroupDialog());
    setupDialogUI(*ui);

    for (auto& [subgroup_type, title] : ::subgrouptype_titles) {
        ui->subgroup_type_combo->addItem(title,
                                         QVariant::fromValue(subgroup_type));
    }

    for (auto& [repeat_pattern, title] : ::repeatpattern_titles) {
        ui->repeat_pattern_combo->addItem(title,
                                          QVariant::fromValue(repeat_pattern));
    }

    m_validator = new NumericLabelValidator(this);

    // Connect signals
    connect(ui->subgroup_type_combo,
            QOverload<int>::of(&QComboBox::currentIndexChanged), this,
            &BracketSubgroupDialog::updateWidgets);

    updateWidgets();
}

AbstractBracketSubgroupDialog::~AbstractBracketSubgroupDialog() = default;

void AbstractBracketSubgroupDialog::updateWidgets()
{
    auto subgroup_type = getSubgroupType();
    bool is_sru = subgroup_type == SubgroupType::SRU_POLYMER;
    QString label_text;
    QString placeholder_text;
    if (is_sru) {
        label_text = "Numeric label:";
        placeholder_text = "(n or value/range: 3, 2-4)";
    } else {
        label_text = "Polymer label:";
    }
    ui->polymer_label_lbl->setText(label_text);

    NumericLabelValidator* validator = is_sru ? m_validator : nullptr;
    ui->polymer_label_le->setValidator(validator);
    ui->polymer_label_le->setPlaceholderText(placeholder_text);
    ui->polymer_label_le->setText(::subgrouptype_label_map.at(subgroup_type));
    ui->polymer_label_le->setEnabled(is_sru);
}

void AbstractBracketSubgroupDialog::setSubgroupType(SubgroupType subgroup_type)
{
    auto idx =
        ui->subgroup_type_combo->findData(QVariant::fromValue(subgroup_type));
    if (idx == -1) {
        // The specified subgroup type is not supported in this combo box
        idx = 0;
    }
    ui->subgroup_type_combo->setCurrentIndex(idx);
}

SubgroupType AbstractBracketSubgroupDialog::getSubgroupType() const
{
    return ui->subgroup_type_combo->currentData().value<SubgroupType>();
}

void AbstractBracketSubgroupDialog::setRepeatPattern(
    RepeatPattern repeat_pattern)
{
    auto idx =
        ui->repeat_pattern_combo->findData(QVariant::fromValue(repeat_pattern));
    if (idx == -1) {
        // The specified repeat pattern is not supported in this combo box
        idx = 0;
    }
    ui->repeat_pattern_combo->setCurrentIndex(idx);
}

RepeatPattern AbstractBracketSubgroupDialog::getRepeatPattern() const
{
    return ui->repeat_pattern_combo->currentData().value<RepeatPattern>();
}

QString AbstractBracketSubgroupDialog::getPolymerLabel() const
{
    return ui->polymer_label_le->text();
}

void AbstractBracketSubgroupDialog::setPolymerLabel(const QString& text)
{
    int pos = 0;
    QString test_text{text};
    auto line_edit = ui->polymer_label_le;
    auto validator = line_edit->validator();
    auto is_valid =
        validator != nullptr &&
        validator->validate(test_text, pos) == QValidator::Acceptable;
    if (line_edit->isEnabled() && is_valid) {
        line_edit->setText(text);
    }
}

BracketSubgroupDialogDeprecated::BracketSubgroupDialogDeprecated(
    QWidget* parent) :
    AbstractBracketSubgroupDialog(parent)
{
}

void BracketSubgroupDialogDeprecated::setAtoms(
    const std::unordered_set<sketcherAtom*>& atoms)
{
    m_atoms = atoms;
}

void BracketSubgroupDialogDeprecated::accept()
{
    emit bracketSubgroupAccepted(getSubgroupType(), getRepeatPattern(),
                                 getPolymerLabel(), m_atoms);

    QDialog::accept();
}

BracketSubgroupDialog::BracketSubgroupDialog(MolModel* const mol_model,
                                             QWidget* parent) :
    AbstractBracketSubgroupDialog(parent),
    m_mol_model(mol_model)
{
}

void BracketSubgroupDialog::accept()
{
    if (auto* atoms = std::get_if<std::unordered_set<const RDKit::Atom*>>(
            &m_atoms_or_s_group)) {
        auto undo_raii = m_mol_model->createUndoMacro("Add substance group");
        m_mol_model->addSGroup(*atoms, getSubgroupType(), getRepeatPattern(),
                               getPolymerLabel().toStdString());
        m_mol_model->clearSelection();
    } else if (auto* s_group = std::get_if<const RDKit::SubstanceGroup*>(
                   &m_atoms_or_s_group)) {
        m_mol_model->modifySGroup(*s_group, getSubgroupType(),
                                  getRepeatPattern(),
                                  getPolymerLabel().toStdString());
    }
    QDialog::accept();
}

void BracketSubgroupDialog::setAtoms(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    m_atoms_or_s_group = atoms;
}

/**
 * @return the relevant data from the given S-group
 * @throw std::out_of_bounds if the S-group type or repeat pattern can not be
 * represented using their respective constants (which means that the settings
 * from this S-group cannot be loaded into this dialog)
 */
static std::tuple<SubgroupType, QString, RepeatPattern>
get_sgroup_data(const RDKit::SubstanceGroup* const s_group)
{
    auto subgroup_type_str = get_sgroup_type(*s_group);
    auto subgroup_type =
        SUBGROUPTYPE_TO_RDKITSTRING_BIMAP.right.at(subgroup_type_str);

    auto polymer_label_str = get_polymer_label(*s_group);
    auto polymer_label = QString::fromStdString(polymer_label_str);

    auto repeat_pattern_str = get_repeat_pattern_label(*s_group);
    auto repeat_pattern =
        REPEATPATTERN_TO_RDKITSTRING_BIMAP.right.at(repeat_pattern_str);

    return {subgroup_type, polymer_label, repeat_pattern};
}

void BracketSubgroupDialog::setSubgroup(
    const RDKit::SubstanceGroup* const s_group)
{
    m_atoms_or_s_group = s_group;

    // load S-group settings into dialog
    SubgroupType subgroup_type;
    QString polymer_label;
    RepeatPattern repeat_pattern;
    try {
        std::tie(subgroup_type, polymer_label, repeat_pattern) =
            get_sgroup_data(s_group);
    } catch (std::out_of_range&) {
        // if the dialog doesn't support this S-group's settings, then leave the
        // dialog settings at their defaults
        return;
    }
    setSubgroupType(subgroup_type);
    setPolymerLabel(polymer_label);
    setRepeatPattern(repeat_pattern);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/bracket_subgroup_dialog.moc"
