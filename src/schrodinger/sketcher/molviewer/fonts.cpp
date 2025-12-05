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
    m_main_label_fm(QFontMetrics(m_main_label_font)),
    m_subscript_fm(QFontMetrics(m_subscript_font)),
    m_charge_fm(QFontMetrics(m_charge_font)),
    m_mapping_fm(QFontMetrics(m_mapping_font)),
    m_chirality_fm(QFontMetrics(m_chirality_font)),
    m_sgroup_fm(QFontMetrics(m_sgroup_font)),
    m_query_label_fm(QFontMetrics(m_query_label_font)),
    m_d_amino_acid_fm(QFontMetrics(m_d_amino_acid_font)),
    m_other_amino_acid_fm(QFontMetrics(m_other_amino_acid_font)),
    m_other_nucleic_acid_sugar_and_phosphate_fm(
        QFontMetrics(m_other_nucleic_acid_sugar_and_phosphate_font)),
    m_other_nucleic_acid_base_fm(QFontMetrics(m_other_nucleic_acid_base_font))
{
    m_d_amino_acid_font.setWeight(QFont::Weight::Bold);
    m_other_amino_acid_font.setWeight(QFont::Weight::Bold);
    m_other_nucleic_acid_sugar_and_phosphate_font.setWeight(
        QFont::Weight::Bold);
    m_other_nucleic_acid_base_font.setWeight(QFont::Weight::Bold);
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
    m_d_amino_acid_font.setPixelSize(size * D_AA_FONT_RATIO);
    m_other_amino_acid_font.setPixelSize(size * OTHER_AA_FONT_RATIO);
    m_other_nucleic_acid_sugar_and_phosphate_font.setPixelSize(
        size * OTHER_NA_SUGAR_AND_PHOSPHATE_FONT_RATIO);
    m_other_nucleic_acid_base_font.setPixelSize(size *
                                                OTHER_NA_BASE_FONT_RATIO);
    m_radical_dot_size = size * RADICAL_DOT_RATIO;
    updateFontMetrics();
}

void Fonts::updateFontMetrics()
{
    // font metrics objects don't update themselves when their font changes, so
    // we need to create new ones
    m_main_label_fm = QFontMetrics(m_main_label_font);
    m_subscript_fm = QFontMetrics(m_subscript_font);
    m_charge_fm = QFontMetrics(m_charge_font);
    m_mapping_fm = QFontMetrics(m_mapping_font);
    m_chirality_fm = QFontMetrics(m_chirality_font);
    m_sgroup_fm = QFontMetrics(m_sgroup_font);
    m_query_label_fm = QFontMetrics(m_query_label_font);

    m_d_amino_acid_fm = QFontMetrics(m_d_amino_acid_font);
    m_other_amino_acid_fm = QFontMetrics(m_other_amino_acid_font);
    m_other_nucleic_acid_sugar_and_phosphate_fm =
        QFontMetrics(m_other_nucleic_acid_sugar_and_phosphate_font);
    m_other_nucleic_acid_base_fm = QFontMetrics(m_other_nucleic_acid_base_font);
}

} // namespace sketcher
} // namespace schrodinger
