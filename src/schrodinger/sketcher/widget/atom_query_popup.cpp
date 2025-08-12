#include "schrodinger/sketcher/widget/atom_query_popup.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_atom_query_popup.h"

namespace schrodinger
{
namespace sketcher
{

AtomQueryPopup::AtomQueryPopup(QWidget* parent) : ModularPopup(parent)
{
    ui.reset(new Ui::AtomQueryPopup());
    ui->setupUi(this);
    setButtonGroup(ui->group);
    setStyleSheet(ATOM_QUERY_STYLE);
}

AtomQueryPopup::~AtomQueryPopup() = default;

void AtomQueryPopup::generateButtonPackets()
{
    m_button_packets.emplace_back(ui->A_btn, static_cast<int>(AtomQuery::A));
    m_button_packets.emplace_back(ui->AH_btn, static_cast<int>(AtomQuery::AH));
    m_button_packets.emplace_back(ui->Q_btn, static_cast<int>(AtomQuery::Q));
    m_button_packets.emplace_back(ui->QH_btn, static_cast<int>(AtomQuery::QH));
    m_button_packets.emplace_back(ui->M_btn, static_cast<int>(AtomQuery::M));
    m_button_packets.emplace_back(ui->MH_btn, static_cast<int>(AtomQuery::MH));
    m_button_packets.emplace_back(ui->X_btn, static_cast<int>(AtomQuery::X));
    m_button_packets.emplace_back(ui->XH_btn, static_cast<int>(AtomQuery::XH));
}

int AtomQueryPopup::getButtonIDToCheck()
{
    auto model = getModel();
    if (model == nullptr) {
        return -1;
    }

    int button_id = -1;
    auto draw_tool = model->getDrawTool();
    auto atom_tool = model->getAtomTool();
    if (draw_tool == DrawTool::ATOM && atom_tool == AtomTool::QUERY) {
        button_id = model->getValueInt(ModelKey::ATOM_QUERY);
    }
    return button_id;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/atom_query_popup.moc"
