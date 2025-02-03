#include "schrodinger/sketcher/molviewer/fonts.h"

#include "schrodinger/sketcher/image_constants.h"
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
    m_query_label_font(QFont(FONT_NAME)),
    m_cursor_hint_font(QFont(FONT_NAME)),
    m_main_label_fm(QFontMetrics(m_main_label_font)),
    m_subscript_fm(QFontMetrics(m_subscript_font)),
    m_charge_fm(QFontMetrics(m_charge_font)),
    m_mapping_fm(QFontMetrics(m_mapping_font)),
    m_chirality_fm(QFontMetrics(m_chirality_font)),
    m_sgroup_fm(QFontMetrics(m_sgroup_font)),
    m_query_label_fm(QFontMetrics(m_query_label_font))
{
    setSize(DEFAULT_FONT_SIZE);
}

int Fonts::size() const
{
    return m_main_label_font.pixelSize();
}

void Fonts::setSize(int size)
{
    if (size == this->size()) {
        return;
    }
    m_main_label_font.setPixelSize(size);
    m_subscript_font.setPixelSize(size * SUBSCRIPT_FONT_RATIO);
    m_charge_font.setPixelSize(size * CHARGE_FONT_RATIO);
    m_mapping_font.setPixelSize(size * MAPPING_FONT_RATIO);
    m_chirality_font.setPixelSize(size * CHIRALITY_FONT_RATIO);
    m_cursor_hint_font.setPixelSize(size * CURSOR_HINT_FONT_RATIO);
    m_sgroup_font.setPixelSize(size * SGROUP_FONT_RATIO);
    m_query_label_font.setPixelSize(size * QUERY_LABEL_FONT_RATIO);

    // font metrics objects don't update themselves when their font changes, so
    // we need to create new ones
    m_main_label_fm = QFontMetrics(m_main_label_font);
    m_subscript_fm = QFontMetrics(m_subscript_font);
    m_charge_fm = QFontMetrics(m_charge_font);
    m_mapping_fm = QFontMetrics(m_mapping_font);
    m_chirality_fm = QFontMetrics(m_chirality_font);
    m_sgroup_fm = QFontMetrics(m_sgroup_font);
    m_query_label_fm = QFontMetrics(m_query_label_font);

    m_radical_dot_size = size * RADICAL_DOT_RATIO;
}

} // namespace sketcher
} // namespace schrodinger
