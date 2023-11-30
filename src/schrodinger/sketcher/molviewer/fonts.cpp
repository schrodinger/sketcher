#include "schrodinger/sketcher/molviewer/fonts.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

Fonts::Fonts() :
    m_main_label_font(QFont(FONT_NAME)),
    m_subscript_font(QFont(FONT_NAME)),
    m_charge_font(QFont(FONT_NAME)),
    m_mapping_font(QFont(FONT_NAME)),
    m_chirality_font(QFont(FONT_NAME)),
    m_sgroup_font(QFont(FONT_NAME)),
    m_cursor_hint_font(QFont(FONT_NAME)),
    m_main_label_fm(QFontMetrics(m_main_label_font)),
    m_subscript_fm(QFontMetrics(m_subscript_font)),
    m_charge_fm(QFontMetrics(m_charge_font)),
    m_mapping_fm(QFontMetrics(m_mapping_font)),
    m_chirality_fm(QFontMetrics(m_chirality_font)),
    m_sgroup_fm(QFontMetrics(m_sgroup_font))
{
    setSize(DEFAULT_FONT_SIZE);
}

qreal Fonts::size() const
{
    return m_main_label_font.pointSizeF();
}

void Fonts::setSize(qreal size)
{
    m_main_label_font.setPointSizeF(size);
    m_subscript_font.setPointSizeF(size * SUBSCRIPT_FONT_RATIO);
    m_charge_font.setPointSizeF(size * CHARGE_FONT_RATIO);
    m_mapping_font.setPointSizeF(size * MAPPING_FONT_RATIO);
    m_chirality_font.setPointSizeF(size * CHIRALITY_FONT_RATIO);
    m_cursor_hint_font.setPointSizeF(size * CURSOR_HINT_FONT_RATIO);
    m_sgroup_font.setPointSizeF(size * SGROUP_FONT_RATIO);

    // font metrics objects don't update themselves when their font changes, so
    // we need to create new ones
    m_main_label_fm = QFontMetrics(m_main_label_font);
    m_subscript_fm = QFontMetrics(m_subscript_font);
    m_charge_fm = QFontMetrics(m_charge_font);
    m_mapping_fm = QFontMetrics(m_mapping_font);
    m_chirality_fm = QFontMetrics(m_chirality_font);
    m_sgroup_fm = QFontMetrics(m_sgroup_font);

    m_radical_dot_size = size * RADICAL_DOT_RATIO;
}

} // namespace sketcher
} // namespace schrodinger
