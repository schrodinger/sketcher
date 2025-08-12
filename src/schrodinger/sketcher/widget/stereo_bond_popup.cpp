#include "schrodinger/sketcher/widget/stereo_bond_popup.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_stereo_bond_popup.h"

namespace schrodinger
{
namespace sketcher
{

StereoBondPopup::StereoBondPopup(QWidget* parent) : AbstractBondPopup(parent)
{
    ui.reset(new Ui::StereoBondPopup());
    ui->setupUi(this);
    setButtonGroup(ui->group);
}

StereoBondPopup::~StereoBondPopup() = default;

void StereoBondPopup::generateButtonPackets()
{
    m_button_packets.emplace_back(ui->up_btn,
                                  static_cast<int>(BondTool::SINGLE_UP));
    m_button_packets.emplace_back(ui->down_btn,
                                  static_cast<int>(BondTool::SINGLE_DOWN));
    m_button_packets.emplace_back(ui->single_either_btn,
                                  static_cast<int>(BondTool::SINGLE_EITHER));
    m_button_packets.emplace_back(ui->double_either_btn,
                                  static_cast<int>(BondTool::DOUBLE_EITHER));
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/stereo_bond_popup.moc"
