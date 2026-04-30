#include "schrodinger/sketcher/widget/nucleic_acid_symbol_popup.h"

#include <QButtonGroup>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/widget/monomer_symbol_popup_utils.h"

namespace schrodinger
{
namespace sketcher
{

NucleicAcidSymbolPopup::NucleicAcidSymbolPopup(
    const std::string& standard_symbol, const std::string& standard_name,
    const std::vector<rdkit_extensions::MonomerInfo>& analogs,
    QWidget* parent) :
    ModularPopup(parent)
{
    auto* group =
        build_monomer_symbol_buttons(this, "na_analog", standard_symbol,
                                     standard_name, analogs, m_id_to_symbol);
    setStyleSheet(ATOM_ELEMENT_OR_MONOMER_STYLE);
    setButtonGroup(group);
}

QString NucleicAcidSymbolPopup::getSymbolForId(int id) const
{
    auto it = m_id_to_symbol.find(id);
    if (it == m_id_to_symbol.end()) {
        return {};
    }
    return QString::fromStdString(it->second);
}

void NucleicAcidSymbolPopup::generateButtonPackets()
{
    auto buttons = m_group->buttons();
    for (int i = 0; i < buttons.size(); ++i) {
        auto* btn = qobject_cast<QToolButton*>(buttons[i]);
        if (btn) {
            m_button_packets.emplace_back(btn, i);
        }
    }
}

int NucleicAcidSymbolPopup::getButtonIDToCheck()
{
    auto model = getModel();
    if (model == nullptr) {
        return -1;
    }

    if (model->getDrawTool() != DrawTool::MONOMER ||
        model->getMonomerToolType() != MonomerToolType::NUCLEIC_ACID) {
        return -1;
    }

    auto analog = model->getValue(ModelKey::NUCLEIC_ACID_SYMBOL)
                      .value<NucleicAcidMutation>();
    for (const auto& [id, symbol] : m_id_to_symbol) {
        if (QString::fromStdString(symbol) == analog.symbol) {
            return id;
        }
    }
    return -1;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/nucleic_acid_symbol_popup.moc"
