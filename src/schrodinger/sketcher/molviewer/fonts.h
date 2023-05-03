#pragma once

#include <QFont>
#include <QFontMetrics>
#include <QtGlobal>

namespace schrodinger
{
namespace sketcher
{

/**
 * An object used for storing all fonts (and associated QFontMetrics) used in a
 * Scene. Also includes the size for radical dots pen.
 */
class Fonts
{
  public:
    Fonts();
    // Prevent implicit copies, since that could cause bugs, as the Scene shares
    // a single Fonts object with all AtomItems.
    explicit Fonts(const Fonts&) = default;
    /**
     * Return the size in points of the font used to paint the main label.  The
     * main label is typically the element symbol, but it can also be a query
     * string or an R group.  The other font sizes are defined as a percentage
     * of the label font size and will be automatically updated.
     */
    qreal size() const;
    void setSize(qreal size);

    QFont m_main_label_font;
    QFont m_subscript_font;
    QFont m_charge_font;
    QFont m_mapping_font;
    QFont m_chirality_font;

    QFontMetrics m_main_label_fm;
    QFontMetrics m_subscript_fm;
    QFontMetrics m_charge_fm;
    QFontMetrics m_mapping_fm;
    QFontMetrics m_chirality_fm;

    qreal m_radical_dot_size;
};

} // namespace sketcher
} // namespace schrodinger
