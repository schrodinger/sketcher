#include "schrodinger/sketcher/molviewer/monomer_label_text.h"

#include <QBrush>
#include <QColor>
#include <QFont>
#include <QFontMetricsF>
#include <QLinearGradient>
#include <QPainter>
#include <QPainterPath>
#include <QPointF>
#include <QRectF>

#include "schrodinger/sketcher/molviewer/monomer_constants.h"

namespace schrodinger
{
namespace sketcher
{

QString elide_text(const std::string& text)
{
    auto qtext = QString::fromStdString(text);
    // A label with one character beyond MAX_MONOMER_LABEL_LENGTH fits without
    // suggesting that more text follows. For longer labels, crop at
    // MAX_MONOMER_LABEL_LENGTH so the additional character can provide the
    // fade-out at the label boundary.
    if (qtext.length() > MAX_MONOMER_LABEL_LENGTH + 1) {
        return qtext.left(MAX_MONOMER_LABEL_LENGTH);
    } else {
        return qtext;
    }
}

QString get_fade_out_text(const QString& text)
{
    if (elide_text(text.toStdString()) == text) {
        return text;
    }
    return text.left(MAX_MONOMER_LABEL_LENGTH + 1);
}

void paint_text_with_fade_out(QPainter& painter, const QString& elided_text,
                              const QString& fade_out_text,
                              const QPointF& left_baseline, const QFont& font,
                              const QColor& color,
                              const QPainterPath& clip_path)
{
    QPainterPath text_path;
    text_path.addText(left_baseline, font, fade_out_text);

    // Start the gradient at the last character used to size the label.
    // Subtracting its advance from the shortened label's advance preserves
    // the kerning Qt applies before that character.
    const QFontMetricsF metrics(font, painter.device());
    const auto fade_start_character = elided_text.right(1);
    const auto fade_start_character_x =
        left_baseline.x() + metrics.horizontalAdvance(elided_text) -
        metrics.horizontalAdvance(fade_start_character);
    QPainterPath fade_start_character_path;
    fade_start_character_path.addText(
        QPointF(fade_start_character_x, left_baseline.y()), font,
        fade_start_character);
    const auto fade_start_x = fade_start_character_path.boundingRect().left();

    // Intersect the glyph outlines explicitly rather than relying only on the
    // painter clip. This ensures that SVG output also truncates the additional
    // character at the clip boundary.
    const auto clipped_text_path = text_path.intersected(clip_path);
    const auto fade_end_x = clipped_text_path.boundingRect().right();

    auto faded_color = color;
    faded_color.setAlpha(0);
    QLinearGradient gradient(fade_start_x, 0.0, fade_end_x, 0.0);
    gradient.setColorAt(0.0, color);
    gradient.setColorAt(1.0, faded_color);
    painter.fillPath(clipped_text_path, QBrush(gradient));
}

QPointF center_text_in_rect(const QString& text, const QFontMetricsF& fm,
                            const QRectF& rect)
{
    auto text_rect = fm.tightBoundingRect(text);
    auto text_rect_origin_loc = -text_rect.bottomLeft();
    text_rect.moveCenter(rect.center());
    return text_rect.bottomLeft() + text_rect_origin_loc;
}

} // namespace sketcher
} // namespace schrodinger
