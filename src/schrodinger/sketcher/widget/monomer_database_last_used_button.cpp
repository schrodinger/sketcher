#include "schrodinger/sketcher/widget/monomer_database_last_used_button.h"

#include <QFontMetricsF>
#include <QPainterPath>
#include <QStyleOptionToolButton>
#include <QStylePainter>

#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/molviewer/monomer_label_text.h"

namespace schrodinger
{
namespace sketcher
{

MonomerDatabaseLastUsedButton::MonomerDatabaseLastUsedButton(QWidget* parent) :
    QToolButton(parent)
{
    setCheckable(true);
    updateMonomerInfo();
}

void MonomerDatabaseLastUsedButton::setMonomer(
    const rdkit_extensions::MonomerID& monomer)
{
    if (monomer == m_monomer) {
        return;
    }
    m_monomer = monomer;
    updateMonomerInfo();
}

const rdkit_extensions::MonomerID&
MonomerDatabaseLastUsedButton::getMonomer() const
{
    return m_monomer;
}

void MonomerDatabaseLastUsedButton::updateMonomerInfo()
{
    const auto symbol = QString::fromStdString(m_monomer.symbol);
    setText(symbol);
    m_elided_text = elide_text(m_monomer.symbol);
    m_fade_out_text = get_fade_out_text(symbol);
    QString tool_tip;
    if (!m_monomer.symbol.empty()) {
        const auto info =
            rdkit_extensions::MonomerDatabase::instance().getMonomerInfo(
                m_monomer.symbol, m_monomer.chain_type);
        tool_tip = QString::fromStdString(info.name.value_or(""));
    }
    setToolTip(tool_tip);
    setStyleSheet(get_atom_element_or_monomer_style(m_elided_text.size()));
    // A blank button looks like empty space, so it mustn't react to the mouse
    setEnabled(!m_monomer.symbol.empty());
    update();
}

bool MonomerDatabaseLastUsedButton::isTextFadedOut() const
{
    return m_elided_text != text();
}

void MonomerDatabaseLastUsedButton::paintEvent(QPaintEvent* event)
{
    if (!isTextFadedOut()) {
        QToolButton::paintEvent(event);
        return;
    }
    // Paint the background without the text, then paint the truncated text
    // centered in the button, the same way that monomer graphics items do
    QStylePainter painter(this);
    QStyleOptionToolButton option;
    initStyleOption(&option);
    option.text.clear();
    painter.drawComplexControl(QStyle::CC_ToolButton, option);

    const QRectF label_rect = rect();
    const auto left_baseline =
        center_text_in_rect(m_elided_text, QFontMetricsF(font()), label_rect);
    QPainterPath clip_path;
    clip_path.addRect(label_rect);
    paint_text_with_fade_out(
        painter, m_elided_text, m_fade_out_text, left_baseline, font(),
        option.palette.color(QPalette::ButtonText), clip_path);
}

void MonomerDatabaseLastUsedButton::nextCheckState()
{
    // The owning tool widget sets the checked state from the model
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/monomer_database_last_used_button.moc"
