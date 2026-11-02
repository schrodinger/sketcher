/**
 * @file
 * @brief Helper functions used by both the C++ tests and
 *        playwright_test_bridge.cpp.
 */

#pragma once

#include <QFontMetricsF>
#include <QGraphicsSimpleTextItem>
#include <QRectF>

namespace schrodinger
{
namespace sketcher
{

/**
 * @return the label's painted text bounds in scene coordinates, excluding the
 * whitespace included in QGraphicsSimpleTextItem::boundingRect().
 */
inline QRectF get_label_ink_rect(const QGraphicsSimpleTextItem& label)
{
    const QFontMetricsF metrics(label.font());
    auto rect = metrics.tightBoundingRect(label.text());
    rect.translate(0.0, metrics.ascent());
    return label.mapRectToScene(rect);
}

} // namespace sketcher
} // namespace schrodinger
