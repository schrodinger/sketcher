#include "schrodinger/sketcher/widget/enumeration_tool_widget.h"

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_enumeration_tool_widget.h"
#include "schrodinger/sketcher/widget/reaction_popup.h"
#include "schrodinger/sketcher/widget/rgroup_popup.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

namespace schrodinger
{
namespace sketcher
{

EnumerationToolWidget::EnumerationToolWidget(QWidget* parent) :
    SketcherView(parent)
{
    ui.reset(new Ui::EnumerationToolWidget());
    ui->setupUi(this);

    m_reaction_popup = new ReactionPopup(this);
    ui->reaction_btn->setPopupWidget(m_reaction_popup);
    ui->reaction_btn->setEnumItem(static_cast<int>(EnumerationTool::RXN_ARROW));

    m_rgroup_popup = new RGroupPopup(this);
    ui->rgroup_btn->setPopupWidget(m_rgroup_popup);
    ui->rgroup_btn->setEnumItem(0);
    ui->rgroup_btn->setStyleSheet(ENUMERATION_STYLE);

    connect(ui->rgroup_btn, &QToolButton::clicked, this,
            &EnumerationToolWidget::onRGroupButtonClicked);
    connect(ui->attachment_point_btn, &QToolButton::clicked, this,
            &EnumerationToolWidget::onAttachmentPointButtonClicked);
    connect(ui->reaction_btn, &QToolButton::clicked, this,
            &EnumerationToolWidget::onReactionButtonClicked);
}

EnumerationToolWidget::~EnumerationToolWidget() = default;

void EnumerationToolWidget::setModel(SketcherModel* model)
{
    SketcherView::setModel(model);
    m_reaction_popup->setModel(model);
    m_rgroup_popup->setModel(model);
}

void EnumerationToolWidget::updateWidgetsEnabled()
{
    auto model = getModel();
    setEnabled(!model->hasActiveSelection());
}

void EnumerationToolWidget::updateCheckState()
{
    updateButtons();

    auto model = getModel();
    auto draw_tool = model->getDrawTool();
    QAbstractButton* button = nullptr;
    if (draw_tool == DrawTool::ENUMERATION) {
        auto enum_tool = model->getEnumerationTool();
        unsigned int rgroup_number =
            model->getValueInt(ModelKey::RGROUP_NUMBER);
        if (enum_tool == EnumerationTool::NEW_RGROUP ||
            (enum_tool == EnumerationTool::EXISTING_RGROUP &&
             m_rgroup_popup->isSupportedRGroup(rgroup_number))) {
            button = ui->rgroup_btn;
        } else if (enum_tool == EnumerationTool::ATTACHMENT_POINT) {
            button = ui->attachment_point_btn;
        } else if (m_reaction_popup->getButtonIDs().count(
                       static_cast<int>(enum_tool)) == 1) {
            button = ui->reaction_btn;
        }
    }

    check_button_or_uncheck_group(button, ui->group);
}

void EnumerationToolWidget::updateButtons()
{
    auto model = getModel();
    auto enum_tool = model->getEnumerationTool();
    if (enum_tool == EnumerationTool::NEW_RGROUP ||
        enum_tool == EnumerationTool::EXISTING_RGROUP) {
        int rgroup_number = enum_tool == EnumerationTool::NEW_RGROUP
                                ? 0u
                                : model->getValueInt(ModelKey::RGROUP_NUMBER);
        ui->rgroup_btn->setEnumItem(rgroup_number);
    } else if (m_reaction_popup->getButtonIDs().count(
                   static_cast<int>(enum_tool)) == 1) {
        ui->reaction_btn->setEnumItem(static_cast<int>(enum_tool));
    }
}

void EnumerationToolWidget::onRGroupButtonClicked()
{
    auto rgroup_number =
        static_cast<unsigned int>(ui->rgroup_btn->getEnumItem());
    auto enum_tool = rgroup_number == 0 ? EnumerationTool::NEW_RGROUP
                                        : EnumerationTool::EXISTING_RGROUP;
    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::ENUMERATION)},
        {ModelKey::ENUMERATION_TOOL, QVariant::fromValue(enum_tool)},
        {ModelKey::RGROUP_NUMBER, QVariant::fromValue(rgroup_number)},
    };
    getModel()->setValues(kv_pairs);
}

void EnumerationToolWidget::onAttachmentPointButtonClicked()
{
    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::ENUMERATION)},
        {ModelKey::ENUMERATION_TOOL,
         QVariant::fromValue(EnumerationTool::ATTACHMENT_POINT)},
    };
    getModel()->setValues(kv_pairs);
}

void EnumerationToolWidget::onReactionButtonClicked()
{
    auto enum_tool = EnumerationTool(ui->reaction_btn->getEnumItem());
    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::ENUMERATION)},
        {ModelKey::ENUMERATION_TOOL, QVariant::fromValue(enum_tool)},
    };
    getModel()->setValues(kv_pairs);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/enumeration_tool_widget.moc"
