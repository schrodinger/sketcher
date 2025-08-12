#include "schrodinger/sketcher/widget/rgroup_popup.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_rgroup_popup.h"

namespace schrodinger
{
namespace sketcher
{

RGroupPopup::RGroupPopup(QWidget* parent) : ModularPopup(parent)
{
    ui.reset(new Ui::RGroupPopup());
    ui->setupUi(this);
    setButtonGroup(ui->group);
    setStyleSheet(ENUMERATION_STYLE);
}

RGroupPopup::~RGroupPopup() = default;

void RGroupPopup::generateButtonPackets()
{
    m_button_packets.emplace_back(ui->new_rgroup_btn, 0);
    m_button_packets.emplace_back(ui->rgroup1_btn, 1);
    m_button_packets.emplace_back(ui->rgroup2_btn, 2);
    m_button_packets.emplace_back(ui->rgroup3_btn, 3);
    m_button_packets.emplace_back(ui->rgroup4_btn, 4);
    m_button_packets.emplace_back(ui->rgroup5_btn, 5);
    m_button_packets.emplace_back(ui->rgroup6_btn, 6);
    m_button_packets.emplace_back(ui->rgroup7_btn, 7);
    m_button_packets.emplace_back(ui->rgroup8_btn, 8);
    m_button_packets.emplace_back(ui->rgroup9_btn, 9);
}

int RGroupPopup::getButtonIDToCheck()
{
    auto model = getModel();
    if (model == nullptr) {
        return -1;
    }

    int button_id = -1;
    auto draw_tool = model->getDrawTool();
    auto enum_tool = model->getEnumerationTool();
    if (draw_tool == DrawTool::ENUMERATION &&
        (enum_tool == EnumerationTool::NEW_RGROUP ||
         enum_tool == EnumerationTool::EXISTING_RGROUP)) {
        auto rgroup_number = model->getValueInt(ModelKey::RGROUP_NUMBER);
        if (enum_tool == EnumerationTool::NEW_RGROUP) {
            button_id = 0;
        } else if (isSupportedRGroup(rgroup_number)) {
            button_id = rgroup_number;
        }
    }
    return button_id;
}

bool RGroupPopup::isSupportedRGroup(unsigned int rgroup_number)
{
    return rgroup_number != 0 && getButtonIDs().count(rgroup_number) == 1;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/rgroup_popup.moc"
