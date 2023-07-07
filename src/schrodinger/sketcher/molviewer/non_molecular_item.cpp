#include "schrodinger/sketcher/molviewer/non_molecular_item.h"

#include <QtGlobal>
#include <QColor>
#include <QPainter>
#include <QPainterPathStroker>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

NonMolecularItem::NonMolecularItem(
    const NonMolecularObject* const non_molecular_object, const QColor& color) :
    m_non_molecular_object(non_molecular_object)
{
    setZValue(static_cast<qreal>(ZOrder::RXN_ARROW_AND_PLUS));
    m_pen = QPen(color);
    m_pen.setWidthF(NON_MOLECULAR_PEN_WIDTH);
    m_pen.setCapStyle(Qt::RoundCap);
    updateCachedData();
}

int NonMolecularItem::type() const
{
    return Type;
}

const NonMolecularObject* NonMolecularItem::getNonMolecularObject() const
{
    return m_non_molecular_object;
}

void NonMolecularItem::updateCachedData()
{
    auto paint_path = QPainterPath();
    // We manually outline the painter path to create the highlight path.  If we
    // used QPainterPathStroker to do this automatically, it would create a path
    // with lots of graphical artifacts from where the path lines cross.
    auto highlight_path = QPainterPath();
    const auto pad = NON_MOLECULAR_HIGHLIGHT_PADDING;
    if (m_non_molecular_object->getType() == NonMolecularType::RXN_ARROW) {
        qreal half_length = ARROW_LENGTH / 2;
        qreal half_width = ARROW_WIDTH / 2;
        qreal tip_start = half_length - half_width;
        paint_path.moveTo(QPointF(-half_length, 0));
        paint_path.lineTo(QPointF(half_length, 0));
        paint_path.lineTo(QPointF(tip_start, -half_width));
        paint_path.lineTo(QPointF(tip_start, half_width));
        paint_path.lineTo(QPointF(half_length, 0));

        // we tweak the padding around the arrow tip so that it looks right
        const auto half_pad = pad / 2.0;
        const auto qtr_pad = pad / 4.0;
        highlight_path.moveTo(QPointF(-half_length - pad, -pad));
        highlight_path.lineTo(QPointF(tip_start - half_pad, -pad));
        highlight_path.lineTo(QPointF(tip_start - half_pad, -pad - half_width));
        highlight_path.lineTo(QPointF(tip_start + qtr_pad, -pad - half_width));
        highlight_path.lineTo(QPointF(half_length + pad + qtr_pad, 0));
        highlight_path.lineTo(QPointF(tip_start + qtr_pad, pad + half_width));
        highlight_path.lineTo(QPointF(tip_start - half_pad, pad + half_width));
        highlight_path.lineTo(QPointF(tip_start - half_pad, pad));
        highlight_path.lineTo(QPointF(-half_length - pad, pad));
        highlight_path.lineTo(QPointF(-half_length - pad, -pad));
    } else { // m_type == NonMolecularType::RXN_PLUS
        qreal arm_length = PLUS_LENGTH / 2;
        paint_path.moveTo(QPointF(-arm_length, 0));
        paint_path.lineTo(QPointF(arm_length, 0));
        paint_path.moveTo(QPointF(0, -arm_length));
        paint_path.lineTo(QPointF(0, arm_length));

        highlight_path.moveTo(QPointF(-arm_length - pad, -pad));
        highlight_path.lineTo(QPointF(-pad, -pad));
        highlight_path.lineTo(QPointF(-pad, -arm_length - pad));
        highlight_path.lineTo(QPointF(pad, -arm_length - pad));
        highlight_path.lineTo(QPointF(pad, -pad));
        highlight_path.lineTo(QPointF(arm_length + pad, -pad));
        highlight_path.lineTo(QPointF(arm_length + pad, pad));
        highlight_path.lineTo(QPointF(pad, pad));
        highlight_path.lineTo(QPointF(pad, arm_length + pad));
        highlight_path.lineTo(QPointF(-pad, arm_length + pad));
        highlight_path.lineTo(QPointF(-pad, pad));
        highlight_path.lineTo(QPointF(-arm_length - pad, pad));
        highlight_path.lineTo(QPointF(-arm_length - pad, -pad));
    }

    m_paint_path = paint_path;
    m_selection_highlighting_path = highlight_path;
    m_predictive_highlighting_path = highlight_path;
    m_shape = highlight_path;
    m_bounding_rect = highlight_path.boundingRect();
}

void NonMolecularItem::paint(QPainter* painter,
                             const QStyleOptionGraphicsItem* option,
                             QWidget* widget)
{
    painter->strokePath(m_paint_path, m_pen);
}

} // namespace sketcher
} // namespace schrodinger
