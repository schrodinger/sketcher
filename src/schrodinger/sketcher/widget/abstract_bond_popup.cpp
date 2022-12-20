#include "schrodinger/sketcher/widget/abstract_bond_popup.h"
#include "schrodinger/sketcher/sketcher_model.h"

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
    auto draw_tool = DrawTool(model->getValue(ModelKey::DRAW_TOOL).toInt());
    if (draw_tool == DrawTool::BOND) {
        button_id = model->getValue(ModelKey::BOND_TOOL).toInt();
    }
    return button_id;
}

} // namespace sketcher
} // namespace schrodinger

#include "abstract_bond_popup.moc"
