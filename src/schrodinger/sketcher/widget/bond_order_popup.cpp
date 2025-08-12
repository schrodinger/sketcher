#include "schrodinger/sketcher/widget/bond_order_popup.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_bond_order_popup.h"

namespace schrodinger
{
namespace sketcher
{

BondOrderPopup::BondOrderPopup(QWidget* parent) : AbstractBondPopup(parent)
{
    ui.reset(new Ui::BondOrderPopup());
    ui->setupUi(this);
    setButtonGroup(ui->group);
}

BondOrderPopup::~BondOrderPopup() = default;

void BondOrderPopup::generateButtonPackets()
{
    m_button_packets.emplace_back(ui->double_btn,
                                  static_cast<int>(BondTool::DOUBLE));
    m_button_packets.emplace_back(ui->triple_btn,
                                  static_cast<int>(BondTool::TRIPLE));
    m_button_packets.emplace_back(ui->coordinate_btn,
                                  static_cast<int>(BondTool::COORDINATE));
    m_button_packets.emplace_back(ui->zero_btn,
                                  static_cast<int>(BondTool::ZERO));
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/bond_order_popup.moc"
