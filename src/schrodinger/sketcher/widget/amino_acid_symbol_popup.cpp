#include "schrodinger/sketcher/widget/amino_acid_symbol_popup.h"

#include <QButtonGroup>
#include <QHBoxLayout>
#include <QToolButton>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"

namespace schrodinger
{
namespace sketcher
{

AminoAcidSymbolPopup::AminoAcidSymbolPopup(
    const std::string& standard_symbol, const std::string& standard_name,
    const std::vector<rdkit_extensions::MonomerInfo>& analogs,
    QWidget* parent) :
    ModularPopup(parent)
{
    auto* layout = new QHBoxLayout(this);
    layout->setContentsMargins(2, 2, 2, 2);
    layout->setSpacing(0);

    auto* group = new QButtonGroup(this);

    // ID 0 = standard amino acid
    int id = 0;
    auto make_button = [&](const std::string& symbol, const std::string& name) {
        auto* btn = new QToolButton(this);
        btn->setText(QString::fromStdString(symbol));
        btn->setToolTip(QString::fromStdString(name));
        btn->setCheckable(true);
        btn->setMinimumSize(32, 30);
        btn->setMaximumSize(32, 32);

        btn->setObjectName(QString::fromStdString("analog_" + symbol + "_btn"));
        btn->setStyleSheet(symbol.size() >= COMPACT_STYLE_MIN_LENGTH
                               ? ATOM_ELEMENT_OR_MONOMER_COMPACT_STYLE
                               : ATOM_ELEMENT_OR_MONOMER_STYLE);
        layout->addWidget(btn);
        group->addButton(btn);
        m_id_to_symbol[id] = symbol;
        ++id;
        return btn;
    };

    make_button(standard_symbol, standard_name);

    for (const auto& analog : analogs) {
        make_button(analog.symbol.value_or(""), analog.name.value_or(""));
    }

    setStyleSheet(ATOM_ELEMENT_OR_MONOMER_STYLE);
    setButtonGroup(group);
}

QString AminoAcidSymbolPopup::getSymbolForId(int id) const
{
    auto it = m_id_to_symbol.find(id);
    if (it == m_id_to_symbol.end()) {
        return {};
    }
    return QString::fromStdString(it->second);
}

void AminoAcidSymbolPopup::generateButtonPackets()
{
    auto buttons = m_group->buttons();
    for (int i = 0; i < buttons.size(); ++i) {
        auto* btn = qobject_cast<QToolButton*>(buttons[i]);
        if (btn) {
            m_button_packets.emplace_back(btn, i);
        }
    }
}

int AminoAcidSymbolPopup::getButtonIDToCheck()
{
    auto model = getModel();
    if (model == nullptr) {
        return -1;
    }

    if (model->getDrawTool() != DrawTool::MONOMER ||
        model->getMonomerToolType() != MonomerToolType::AMINO_ACID) {
        return -1;
    }

    auto analog = model->getValueString(ModelKey::AMINO_ACID_SYMBOL);
    for (const auto& [id, symbol] : m_id_to_symbol) {
        if (QString::fromStdString(symbol) == analog) {
            return id;
        }
    }
    return -1;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/amino_acid_symbol_popup.moc"
