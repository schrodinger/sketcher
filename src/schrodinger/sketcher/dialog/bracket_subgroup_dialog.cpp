#include "schrodinger/sketcher/dialog/bracket_subgroup_dialog.h"

#include <iostream>
#include <unordered_map>

#include <QString>

#include "schrodinger/sketcher/Atom.h"
#include "schrodinger/sketcher/numeric_label_validator.h"
#include "schrodinger/sketcher/s_group_constants.h"
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

BracketSubgroupDialog::BracketSubgroupDialog(QWidget* parent) :
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

BracketSubgroupDialog::~BracketSubgroupDialog() = default;

void BracketSubgroupDialog::updateWidgets()
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

void BracketSubgroupDialog::setSubgroupType(SubgroupType subgroup_type)
{
    auto idx =
        ui->subgroup_type_combo->findData(QVariant::fromValue(subgroup_type));
    if (idx == -1) {
        // The specified subgroup type is not supported in this combo box
        idx = 0;
    }
    ui->subgroup_type_combo->setCurrentIndex(idx);
}

SubgroupType BracketSubgroupDialog::getSubgroupType() const
{
    return ui->subgroup_type_combo->currentData().value<SubgroupType>();
}

void BracketSubgroupDialog::setRepeatPattern(RepeatPattern repeat_pattern)
{
    auto idx =
        ui->repeat_pattern_combo->findData(QVariant::fromValue(repeat_pattern));
    if (idx == -1) {
        // The specified repeat pattern is not supported in this combo box
        idx = 0;
    }
    ui->repeat_pattern_combo->setCurrentIndex(idx);
}

RepeatPattern BracketSubgroupDialog::getRepeatPattern() const
{
    return ui->repeat_pattern_combo->currentData().value<RepeatPattern>();
}

void BracketSubgroupDialog::accept()
{
    emit bracketSubgroupAccepted(getSubgroupType(), getRepeatPattern(),
                                 getPolymerLabel(), m_atoms);

    QDialog::accept();
}

QString BracketSubgroupDialog::getPolymerLabel() const
{
    return ui->polymer_label_le->text();
}

void BracketSubgroupDialog::setPolymerLabel(const QString& text)
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

void BracketSubgroupDialog::setAtoms(std::unordered_set<sketcherAtom*> atoms)
{
    m_atoms = atoms;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/bracket_subgroup_dialog.moc"
