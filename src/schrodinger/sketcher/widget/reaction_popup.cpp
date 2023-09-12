#include "schrodinger/sketcher/widget/reaction_popup.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_reaction_popup.h"

namespace schrodinger
{
namespace sketcher
{

ReactionPopup::ReactionPopup(QWidget* parent) : ModularPopup(parent)
{
    ui.reset(new Ui::ReactionPopup());
    ui->setupUi(this);
    setButtonGroup(ui->group);
}

ReactionPopup::~ReactionPopup() = default;

void ReactionPopup::generateButtonPackets()
{
    m_button_packets.emplace_back(ui->rxn_arrow_btn,
                                  static_cast<int>(EnumerationTool::RXN_ARROW));
    m_button_packets.emplace_back(ui->rxn_plus_btn,
                                  static_cast<int>(EnumerationTool::RXN_PLUS));
    m_button_packets.emplace_back(
        ui->add_mapping_btn, static_cast<int>(EnumerationTool::ADD_MAPPING));
    m_button_packets.emplace_back(
        ui->remove_mapping_btn,
        static_cast<int>(EnumerationTool::REMOVE_MAPPING));
}

int ReactionPopup::getButtonIDToCheck()
{
    auto model = getModel();
    if (model == nullptr) {
        return -1;
    }

    int button_id = -1;
    auto draw_tool = model->getDrawTool();
    if (draw_tool == DrawTool::ENUMERATION) {
        auto enum_tool_int = model->getValueInt(ModelKey::ENUMERATION_TOOL);
        button_id = getButtonIDs().count(enum_tool_int) == 1 ? enum_tool_int
                                                             : button_id;
    }
    return button_id;
}

void ReactionPopup::updateWidgetsEnabled()
{
    ModularPopup::updateWidgetsEnabled();
    auto model = getModel();
    auto allow_multiple_rxns =
        model->getValueBool(ModelKey::ALLOW_MULTIPLE_RXNS);
    ui->rxn_arrow_btn->setEnabled(allow_multiple_rxns || !model->hasReaction());
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/reaction_popup.moc"
