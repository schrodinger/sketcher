#include "schrodinger/sketcher/widget/modular_popup.h"

#include <stdexcept>

#include <QButtonGroup>
#include <QIcon>
#include <QPainter>
#include <QStyleOption>
#include <QVariant>

#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

namespace schrodinger
{
namespace sketcher
{

ButtonPacket::ButtonPacket(QToolButton* button, int enum_int) :
    button(button),
    enum_int(enum_int)
{
}

ModularPopup::ModularPopup(QWidget* parent) :
    SketcherViewWithWasmOutlineFix(parent)
{
    // Show as a pop-up top-level window
    setWindowFlags(Qt::Popup);
}

void ModularPopup::setButtonGroup(QButtonGroup* group)
{
    if (m_group != nullptr) {
        throw std::runtime_error(
            "The button group has already been set on this popup.");
    }

    m_group = group;
    generateButtonPackets();

    // Assign each button its ID based on the specified enum items
    for (auto& packet : m_button_packets) {
        group->setId(packet.button, packet.enum_int);
    }

    connect(group, &QButtonGroup::idClicked, this,
            &ModularPopup::onButtonClicked);
}

QIcon ModularPopup::getIcon(int enum_int) const
{
    auto button = m_group->button(enum_int);
    if (button == nullptr) {
        // No button for the specified enum item, so return a blank icon
        return QIcon();
    }
    return button->icon();
}

QString ModularPopup::getText(int enum_int) const
{
    auto button = m_group->button(enum_int);
    return button == nullptr ? QString() : button->text();
}

QString ModularPopup::getToolTip(int enum_int) const
{
    auto button = m_group->button(enum_int);
    return button == nullptr ? QString()
                             : button->toolTip() + " – press & hold to change";
}

std::vector<ButtonPacket> ModularPopup::getButtonPackets() const
{
    return m_button_packets;
}

void ModularPopup::onButtonClicked(int button_id)
{
    emit selectionChanged(button_id);
    close();
}

void ModularPopup::updateCheckState()
{
    if (m_group == nullptr) {
        return;
    }
    int button_id = (getModel() == nullptr) ? -1 : getButtonIDToCheck();
    auto button = m_group->button(button_id);
    check_button_or_uncheck_group(button, m_group);
}

std::unordered_set<int> ModularPopup::getButtonIDs() const
{
    std::unordered_set<int> button_ids;
    if (m_group != nullptr) {
        for (auto button : m_group->buttons()) {
            button_ids.insert(m_group->id(button));
        }
    }
    return button_ids;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/modular_popup.moc"
