#include "schrodinger/sketcher/widget/abstract_bond_popup.h"

#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

AbstractBondPopup::AbstractBondPopup(QWidget* parent) : ModularPopup(parent)
{
}

int AbstractBondPopup::getButtonIDToCheck()
{
    auto model = getModel();
    if (model == nullptr) {
        return -1;
    }

    int button_id = -1;
    auto draw_tool = model->getDrawTool();
    if (draw_tool == DrawTool::BOND) {
        button_id = model->getValueInt(ModelKey::BOND_TOOL);
    }
    return button_id;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/abstract_bond_popup.moc"
