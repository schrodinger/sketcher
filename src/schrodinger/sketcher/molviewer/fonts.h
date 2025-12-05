#pragma once

#include <QFont>
#include <QFontMetricsF>
#include <QtGlobal>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * An object used for storing all fonts (and associated QFontMetrics) used in a
 * Scene. Also includes the size for radical dots pen.
 */
class SKETCHER_API Fonts
{
  public:
    Fonts();
    // Prevent implicit copies, since that could cause bugs, as the Scene shares
    // a single Fonts object with all AtomItems.
    explicit Fonts(const Fonts&) = default;

    /**
     * Return the size in pixels of the font used to paint the main label.  The
     * main label is typically the element symbol, but it can also be a query
     * string or an R group.  The other font sizes are defined as a percentage
     * of the label font size and will be automatically updated.
     */
    int size() const;
    void setSize(int size);

    /**
     * Update the font metrics after manually making changes to this object's
     * fonts.
     *
     * @note This method is called automatically by setSize() and only needs to
     * be explicitly called when directly manipulating the individual font
     * objects.
     */
    void updateFontMetrics();

    QFont m_main_label_font;
    QFont m_subscript_font;
    QFont m_charge_font;
    QFont m_mapping_font;
    QFont m_chirality_font;
    QFont m_sgroup_font;
    QFont m_query_label_font;
    QFont m_cursor_hint_font;

    QFont m_d_amino_acid_font;
    QFont m_other_amino_acid_font;
    QFont m_other_nucleic_acid_sugar_and_phosphate_font;
    QFont m_other_nucleic_acid_base_font;

    QFontMetricsF m_main_label_fm;
    QFontMetricsF m_subscript_fm;
    QFontMetricsF m_charge_fm;
    QFontMetricsF m_mapping_fm;
    QFontMetricsF m_chirality_fm;
    QFontMetricsF m_sgroup_fm;
    QFontMetricsF m_query_label_fm;
    // we don't need a QFontMetricsF for the cursor hint, so it's omitted here

    QFontMetricsF m_d_amino_acid_fm;
    QFontMetricsF m_other_amino_acid_fm;
    QFontMetricsF m_other_nucleic_acid_sugar_and_phosphate_fm;
    QFontMetricsF m_other_nucleic_acid_base_fm;

    qreal m_radical_dot_size;
};

} // namespace sketcher
} // namespace schrodinger
