#pragma once

#include <string>

#include <QString>

#include "schrodinger/sketcher/definitions.h"

class QColor;
class QFont;
class QPainter;
class QPainterPath;
class QFontMetricsF;
class QPointF;
class QRectF;

namespace schrodinger
{
namespace sketcher
{

/**
 * Return labels containing at most one character beyond
 * MAX_MONOMER_LABEL_LENGTH unchanged. Crop longer labels to
 * MAX_MONOMER_LABEL_LENGTH so one additional character can be painted as the
 * faded truncation character (see paint_text_with_fade_out()).
 */
SKETCHER_API QString elide_text(const std::string& text);

/**
 * @return The text to paint for a label that elide_text() truncated: the
 * truncated label plus the faded truncation character. Returns the label
 * unchanged if elide_text() wouldn't truncate it.
 */
SKETCHER_API QString get_fade_out_text(const QString& text);

/**
 * Paint a label that elide_text() truncated, fading out the label from its
 * last displayed character through the truncation character. The text is
 * clipped to the given path.
 *
 * Note that the label should be sized and centered using elided_text, so the
 * truncation character extends into the label's right-hand padding.
 *
 * @param painter The painter to use
 * @param elided_text The truncated label, as returned by elide_text()
 * @param fade_out_text The text to paint, as returned by get_fade_out_text()
 * @param left_baseline The leftmost point of the baseline to paint the text at
 * @param font The font to paint the text with
 * @param color The text color before it starts fading out
 * @param clip_path The text is clipped to this path, e.g. the label's border
 */
SKETCHER_API void
paint_text_with_fade_out(QPainter& painter, const QString& elided_text,
                         const QString& fade_out_text,
                         const QPointF& left_baseline, const QFont& font,
                         const QColor& color, const QPainterPath& clip_path);

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
SKETCHER_API QPointF center_text_in_rect(const QString& text,
                                         const QFontMetricsF& fm,
                                         const QRectF& rect);

} // namespace sketcher
} // namespace schrodinger
