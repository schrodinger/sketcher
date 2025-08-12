#include "schrodinger/sketcher/widget/modular_tool_button.h"

#include <stdexcept>

#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace schrodinger
{
namespace sketcher
{

ModularToolButton::ModularToolButton(QWidget* parent) :
    ToolButtonWithPopup(parent)
{
}

void ModularToolButton::setPopupWidget(QWidget* popup)
{
    ToolButtonWithPopup::setPopupWidget(popup);
    auto cast_popup = dynamic_cast<ModularPopup*>(m_popup_wdg);
    if (cast_popup == nullptr) {
        throw std::runtime_error("The popup widget for a modular tool button"
                                 " must be a ModularPopup.");
    }
    connect(cast_popup, &ModularPopup::selectionChanged, this,
            &ModularToolButton::onPopupSelectionChanged);
    updateButton();
}

int ModularToolButton::getEnumItem()
{
    return m_enum_int;
}

void ModularToolButton::onPopupSelectionChanged(int enum_int)
{
    setEnumItem(enum_int);
    click();
}

void ModularToolButton::setEnumItem(int enum_int)
{
    if (enum_int != m_enum_int) {
        m_enum_int = enum_int;
        updateButton();
        emit enumChanged();
    }
}

void ModularToolButton::updateButton()
{
    if (m_popup_wdg == nullptr) {
        return;
    }
    auto popup = dynamic_cast<ModularPopup*>(m_popup_wdg);
    setIcon(popup->getIcon(m_enum_int));
    setText(popup->getText(m_enum_int));
    setToolTip(popup->getToolTip(m_enum_int));
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/modular_tool_button.moc"
