#include "schrodinger/sketcher/widget/bond_query_popup.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_bond_query_popup.h"

namespace schrodinger
{
namespace sketcher
{

BondQueryPopup::BondQueryPopup(QWidget* parent) : AbstractBondPopup(parent)
{
    ui.reset(new Ui::BondQueryPopup());
    ui->setupUi(this);
    setButtonGroup(ui->group);
    setStyleSheet(BOND_QUERY_STYLE);
}

BondQueryPopup::~BondQueryPopup() = default;

void BondQueryPopup::generateButtonPackets()
{
    m_button_packets.emplace_back(ui->aromatic_btn,
                                  static_cast<int>(BondTool::AROMATIC));
    m_button_packets.emplace_back(ui->any_btn, static_cast<int>(BondTool::ANY));
    m_button_packets.emplace_back(ui->single_double_btn,
                                  static_cast<int>(BondTool::SINGLE_OR_DOUBLE));
    m_button_packets.emplace_back(
        ui->single_aromatic_btn,
        static_cast<int>(BondTool::SINGLE_OR_AROMATIC));
    m_button_packets.emplace_back(
        ui->double_aromatic_btn,
        static_cast<int>(BondTool::DOUBLE_OR_AROMATIC));
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/bond_query_popup.moc"
