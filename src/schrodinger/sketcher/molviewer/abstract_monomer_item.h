#pragma once

#include <string>
#include <tuple>

#include <QBrush>
#include <QFont>
#include <QPainterPath>
#include <QPen>
#include <QRectF>
#include <QString>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

namespace RDKit
{
class Atom;
}

namespace schrodinger
{
namespace rdkit_extensions
{
enum class ChainType;
}

namespace sketcher
{

/**
 * An abstract parent class for graphics items that represent monomers.
 */
class SKETCHER_API AbstractMonomerItem : public AbstractAtomOrMonomerItem
{
  public:
    AbstractMonomerItem(const RDKit::Atom* monomer, const Fonts& fonts,
                        QGraphicsItem* parent = nullptr);

    /**
     * Replace the monomer-specific colors with the given colors.
     */
    void setMonomerColors(const QColor& background_color,
                          const QColor& outline_color,
                          const QColor& font_color);

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

  protected:
    const Fonts& m_fonts;
    QPen m_border_pen;
    QBrush m_border_brush = QBrush(Qt::BrushStyle::SolidPattern);
    QString m_main_label_text;
    QFont m_main_label_font;
    QPen m_main_label_pen =
        QPen(MONOMER_LABEL_TEXT_COLOR, MONOMER_LABEL_TEXT_WIDTH);
    QPainterPath m_border_path;
    // The leftmost point of the baseline to use when painting the main label
    // text (i.e. what we should pass to the QPainter in order to center the
    // main label)
    QPointF m_main_label_left_baseline;

    /**
     * @return the number scaled the same as the current font setting. For
     * example, if the user has increased the font size by 50% (e.g. 27 pt
     * instead of the default 18 pt), then return num * 1.5.
     */
    qreal scaleBasedOnFontSize(qreal num) const;
};

/**
 * Elide the given string so that it contains no more than 6 characters.
 */
SKETCHER_API QString elide_text(const std::string& text);

/**
 * Return the residue name of the given monomer
 * @param monomer an atom that represents a monomer
 */
SKETCHER_API std::string get_monomer_res_name(const RDKit::Atom* const monomer);

/**
 * Return whether the given residue name is considered a standard nucleotide
 * (i.e. A, T, G, C, or U)
 */
SKETCHER_API bool is_standard_nucleotide(const std::string& res_name);

/**
 * Determine the size of rectangle required to fit the given label text
 * @param main_label_text the text to size the rectangle for
 * @param fm a font metrics object for the font that will be used to display the
 * text
 * @param standard_width the minimum width of the rectangle
 * @param standard_height the minimum height of the rectangle
 */
SKETCHER_API std::pair<qreal, qreal>
get_rect_size_to_fit_label(const QString& main_label_text,
                           const QFontMetricsF& fm, const qreal standard_width,
                           const qreal standard_height);

/**
 * Determine the size of diamond required to fit the given label text
 * @param main_label_text the text to size the diamond for
 * @param fm a font metrics object for the font that will be used to display the
 * text
 * @param standard_width the minimum width of the diamond
 * @param standard_height the minimum height of the diamond
 */
SKETCHER_API std::pair<qreal, qreal> get_diamond_size_to_fit_label(
    const QString& main_label_text, const QFontMetricsF& fm,
    const qreal standard_width, const qreal standard_height);

/**
 * Determine the size of ellipse required to fit the given label text
 * @param main_label_text the text to size the ellipse for
 * @param fm a font metrics object for the font that will be used to display the
 * text
 * @param standard_width the minimum width of the ellipse
 * @param standard_height the minimum height of the ellipse
 */
SKETCHER_API std::pair<qreal, qreal> get_ellipse_size_to_fit_label(
    const QString& main_label_text, const QFontMetricsF& fm,
    const qreal standard_width, const qreal standard_height);

/**
 * Get the border pen thickness, border color, and font to use (along with a
 * font metrics object) for drawing a nucleic acid monomer of the given
 * residue name.
 * @param res_name the name of the residue being drawn
 * @param fonts the fonts object to pull fonts and font metrics objects from
 * @param use_base_font if true, return a font for nucleobases.  Otherwise,
 * return a font for sugars and phosphates
 */
SKETCHER_API std::tuple<qreal, QColor, const QFont, const QFontMetricsF*>
get_border_and_font_settings_for_nucleic_acid(const std::string& res_name,
                                              const Fonts& fonts,
                                              const bool use_base_font = false);

/**
 * Update the given path so it contains only a rectangle of the specified size
 * @param[in,out] path the path to modify
 * @param[in] rect the rectangle to add to the path
 * @param[in] highlighting_thickness additional size to add to the rectangle.
 * This option is intended for generating selection and predictive highlighting
 * paths.
 */
SKETCHER_API void set_path_to_rect(QPainterPath& path, const QRectF& rect,
                                   const qreal highlighting_thickness = 0);

/**
 * @return a new rect that's larger than the given rect by half of the specified
 * pen width on all four sides
 *
 * @note This is typically used to create a bounding rect that accounts for the
 * width of the pen used to draw an outline. The drawn line is centered on
 * the outline, so half of the pen width is *outside* of the outline.
 */
QRectF rect_expanded_by_half_pen_width(const QRectF& rect,
                                       const qreal pen_width);

/**
 * Get the appropriate color to use for a monomer of the given residue name.
 * If the residue name is not in the color map, this function will check the
 * monomer database for a natural analog and use its color if available.
 * @param res_name the residue name of the monomer
 * @param chain_type the chain type of the monomer (PEPTIDE, RNA, DNA, or CHEM)
 * @param color_by_res_name a mapping of residue name to color
 * @param default_color the default color to use if the residue name is not
 * found in the mapping and no natural analog is available
 */
SKETCHER_API QColor get_color_for_monomer(
    const std::string& res_name, rdkit_extensions::ChainType chain_type,
    const std::unordered_map<std::string, QColor>& color_by_res_name,
    const QColor& default_color);

/**
 * Determine the coordinates of the left end of the base line required to center
 * the specified text in the given rectangle. In other words, figure out what
 * point to pass to QPainter::drawText(QPointF, QString) in order to center the
 * text.
 * @param text The text to center
 * @param fm The font metrics for the font that will be used to render the text
 * @param rect The rectangle to center the text in
 * @return The base line coordinate
 * @note If we instead painted our text by calling QPainter::drawText(rect,
 * Qt::AlignCenter, text), then the vertical centering would be based on the
 * height of a "typical" letter (i.e. cap height / 2), *not* the height of the
 * actual text to be rendered. This is particularly noticeable with the capital
 * Q for a glutamine monomer, where the descender (i.e. the bottom bit of the
 * letter that sits below the baseline) is completely ignored for centering
 * purposes. With large font sizes, this leads to the descender sticking out of
 * the monomer's shape even though there's a fair bit of empty space above the
 * letter. This function avoids that problem by centering the text using the
 * actual size of the letters.
 */
QPointF center_text_in_rect(const QString& text, const QFontMetricsF& fm,
                            const QRectF& rect);

} // namespace sketcher
} // namespace schrodinger
