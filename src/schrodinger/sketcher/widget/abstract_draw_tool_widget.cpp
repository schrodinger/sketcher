#include "schrodinger/sketcher/widget/abstract_draw_tool_widget.h"

#include <QAbstractButton>

#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

AbstractDrawToolWidget::AbstractDrawToolWidget(QWidget* parent) :
    SketcherView(parent)
{
}

std::unordered_set<QAbstractButton*>
AbstractDrawToolWidget::getCheckableButtons()
{
    return std::unordered_set<QAbstractButton*>();
}

void AbstractDrawToolWidget::updateCheckState()
{
    auto model = getModel();
    if (model == nullptr) {
        return;
    }
    bool checkable = !model->hasActiveSelection();
    for (auto button : getCheckableButtons()) {
        button->setCheckable(checkable);
    }
    updateCheckedButton();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/abstract_draw_tool_widget.moc"
