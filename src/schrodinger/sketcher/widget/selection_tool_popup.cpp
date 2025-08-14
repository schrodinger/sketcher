#include "schrodinger/sketcher/ui/ui_selection_tool_popup.h"
#include "schrodinger/sketcher/widget/selection_tool_popup.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

SelectionToolPopup::SelectionToolPopup(QWidget* parent) : ModularPopup(parent)
{
    ui.reset(new Ui::SelectionToolPopup());
    ui->setupUi(this);
    setButtonGroup(ui->group);
}

SelectionToolPopup::~SelectionToolPopup()
{
}

int SelectionToolPopup::getButtonIDToCheck()
{
    auto model = getModel();
    if (model == nullptr) {
        return -1;
    }

    int button_id = -1;
    auto draw_tool = model->getDrawTool();
    if (draw_tool == DrawTool::SELECT) {
        button_id = model->getValueInt(ModelKey::SELECTION_TOOL);
    }
    return button_id;
}

void SelectionToolPopup::generateButtonPackets()
{
    m_button_packets.emplace_back(ui->rect_btn,
                                  static_cast<int>(SelectionTool::RECTANGLE));
    m_button_packets.emplace_back(ui->lasso_btn,
                                  static_cast<int>(SelectionTool::LASSO));
    m_button_packets.emplace_back(ui->ellipse_btn,
                                  static_cast<int>(SelectionTool::ELLIPSE));
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/selection_tool_popup.moc"
