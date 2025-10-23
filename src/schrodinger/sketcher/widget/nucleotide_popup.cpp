#include "schrodinger/sketcher/ui/ui_nucleotide_popup.h"

#include <QButtonGroup>

#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/widget/nucleotide_popup.h"

namespace schrodinger
{
namespace sketcher
{

NucleotidePopup::NucleotidePopup(const NucleicAcidTool tool,
                                 const ModelKey model_key, const QString& sugar,
                                 const QString& u_or_t, QWidget* parent) :
    ModularPopup(parent),
    m_tool(tool),
    m_model_key(model_key),
    m_sugar(sugar),
    m_u_or_t(u_or_t)
{
    ui.reset(new Ui::NucleotidePopup());
    ui->setupUi(this);
    setButtonGroup(ui->group);
    setStyleSheet(ATOM_ELEMENT_OR_MONOMER_STYLE);

    // add text to the buttons
    QString btn_name_fmt("%1(%2)P");
    ui->a_btn->setText(btn_name_fmt.arg(sugar, "A"));
    ui->c_btn->setText(btn_name_fmt.arg(sugar, "C"));
    ui->g_btn->setText(btn_name_fmt.arg(sugar, "G"));
    ui->u_or_t_btn->setText(btn_name_fmt.arg(sugar, u_or_t));
    ui->n_btn->setText(btn_name_fmt.arg(sugar, "N"));
}

NucleotidePopup::~NucleotidePopup() = default;

int NucleotidePopup::getButtonIDToCheck()
{
    auto model = getModel();
    if (model == nullptr) {
        return -1;
    }

    int button_id = -1;
    if (model->getDrawTool() == DrawTool::MONOMER &&
        model->getMonomerToolType() == MonomerToolType::NUCLEIC_ACID &&
        model->getNucleicAcidTool() == m_tool) {
        StdNucleobase base =
            model->getValue(m_model_key).value<StdNucleobase>();
        switch (base) {
            case StdNucleobase::A:
                button_id = m_group->id(ui->a_btn);
                break;
            case StdNucleobase::C:
                button_id = m_group->id(ui->c_btn);
                break;
            case StdNucleobase::G:
                button_id = m_group->id(ui->g_btn);
                break;
            case StdNucleobase::U_OR_T:
                button_id = m_group->id(ui->u_or_t_btn);
                break;
            case StdNucleobase::N:
                button_id = m_group->id(ui->n_btn);
                break;
            default:
                break;
        }
    }
    return button_id;
}

void NucleotidePopup::generateButtonPackets()
{
    m_button_packets.emplace_back(ui->a_btn,
                                  static_cast<int>(StdNucleobase::A));
    m_button_packets.emplace_back(ui->c_btn,
                                  static_cast<int>(StdNucleobase::C));
    m_button_packets.emplace_back(ui->g_btn,
                                  static_cast<int>(StdNucleobase::G));
    m_button_packets.emplace_back(ui->u_or_t_btn,
                                  static_cast<int>(StdNucleobase::U_OR_T));
    m_button_packets.emplace_back(ui->n_btn,
                                  static_cast<int>(StdNucleobase::N));
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/nucleotide_popup.moc"
