#include "schrodinger/sketcher/molviewer/fonts.h"

#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"

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
    m_d_amino_acid_font(QFont(FONT_NAME)),
    m_other_amino_acid_font(QFont(FONT_NAME)),
    m_other_nucleic_acid_sugar_and_phosphate_font(QFont(FONT_NAME)),
    m_other_nucleic_acid_base_font(QFont(FONT_NAME)),
    m_monomeric_attachment_point_label_font(QFont(FONT_NAME)),
    m_main_label_fm(QFontMetricsF(m_main_label_font)),
    m_subscript_fm(QFontMetricsF(m_subscript_font)),
    m_charge_fm(QFontMetricsF(m_charge_font)),
    m_mapping_fm(QFontMetricsF(m_mapping_font)),
    m_chirality_fm(QFontMetricsF(m_chirality_font)),
    m_sgroup_fm(QFontMetricsF(m_sgroup_font)),
    m_query_label_fm(QFontMetricsF(m_query_label_font)),
    m_d_amino_acid_fm(QFontMetricsF(m_d_amino_acid_font)),
    m_other_amino_acid_fm(QFontMetricsF(m_other_amino_acid_font)),
    m_other_nucleic_acid_sugar_and_phosphate_fm(
        QFontMetricsF(m_other_nucleic_acid_sugar_and_phosphate_font)),
    m_other_nucleic_acid_base_fm(QFontMetricsF(m_other_nucleic_acid_base_font)),
    m_monomeric_attachment_point_label_fm(
        QFontMetricsF(m_monomeric_attachment_point_label_font))
{
    m_d_amino_acid_font.setWeight(QFont::Weight::Bold);
    m_other_amino_acid_font.setWeight(QFont::Weight::Bold);
    m_other_nucleic_acid_sugar_and_phosphate_font.setWeight(
        QFont::Weight::Bold);
    m_other_nucleic_acid_base_font.setWeight(QFont::Weight::Bold);
    m_monomeric_attachment_point_label_font.setWeight(QFont::Weight::Bold);
    setSize(DEFAULT_FONT_SIZE);
}

int Fonts::size() const
{
    return m_main_label_font.pixelSize();
}

/**
 * Call setPixelSize on the font after ensuring that the size is at least 1
 * pixel (since calling setPixelSize(0) is a no-op and prints a warning to the
 * terminal)
 */
static void set_font_pixel_size(QFont& font, int size)
{
    size = std::max(size, 1);
    font.setPixelSize(size);
}

void Fonts::setSize(int size)
{
    if (size <= 0 || size == this->size()) {
        return;
    }
    set_font_pixel_size(m_main_label_font, size);
    set_font_pixel_size(m_subscript_font, size * SUBSCRIPT_FONT_RATIO);
    set_font_pixel_size(m_charge_font, size * CHARGE_FONT_RATIO);
    set_font_pixel_size(m_mapping_font, size * MAPPING_FONT_RATIO);
    set_font_pixel_size(m_chirality_font, size * CHIRALITY_FONT_RATIO);
    set_font_pixel_size(m_cursor_hint_font, size * CURSOR_HINT_FONT_RATIO);
    set_font_pixel_size(m_sgroup_font, size * SGROUP_FONT_RATIO);
    set_font_pixel_size(m_query_label_font, size * QUERY_LABEL_FONT_RATIO);
    set_font_pixel_size(m_d_amino_acid_font, size * D_AA_FONT_RATIO);
    set_font_pixel_size(m_other_amino_acid_font, size * OTHER_AA_FONT_RATIO);
    set_font_pixel_size(m_other_nucleic_acid_sugar_and_phosphate_font,
                        size * OTHER_NA_SUGAR_AND_PHOSPHATE_FONT_RATIO);
    set_font_pixel_size(m_other_nucleic_acid_base_font,
                        size * OTHER_NA_BASE_FONT_RATIO);
    set_font_pixel_size(m_monomeric_attachment_point_label_font,
                        size * MONOMERIC_ATTACHMENT_POINT_LABEL_FONT_RATIO);
    m_radical_dot_size = size * RADICAL_DOT_RATIO;
    updateFontMetrics();
}

void Fonts::updateFontMetrics()
{
    // font metrics objects don't update themselves when their font changes, so
    // we need to create new ones
    m_main_label_fm = QFontMetricsF(m_main_label_font);
    m_subscript_fm = QFontMetricsF(m_subscript_font);
    m_charge_fm = QFontMetricsF(m_charge_font);
    m_mapping_fm = QFontMetricsF(m_mapping_font);
    m_chirality_fm = QFontMetricsF(m_chirality_font);
    m_sgroup_fm = QFontMetricsF(m_sgroup_font);
    m_query_label_fm = QFontMetricsF(m_query_label_font);

    m_d_amino_acid_fm = QFontMetricsF(m_d_amino_acid_font);
    m_other_amino_acid_fm = QFontMetricsF(m_other_amino_acid_font);
    m_other_nucleic_acid_sugar_and_phosphate_fm =
        QFontMetricsF(m_other_nucleic_acid_sugar_and_phosphate_font);
    m_other_nucleic_acid_base_fm =
        QFontMetricsF(m_other_nucleic_acid_base_font);
    m_monomeric_attachment_point_label_fm =
        QFontMetricsF(m_monomeric_attachment_point_label_font);
}

} // namespace sketcher
} // namespace schrodinger
