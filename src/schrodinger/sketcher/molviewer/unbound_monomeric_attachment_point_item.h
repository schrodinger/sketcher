#pragma once

#include <QBrush>
#include <QGraphicsItem>
#include <QPen>
#include <QPointF>
#include <QRectF>
#include <QString>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace schrodinger
{
namespace sketcher
{

class AbstractMonomerItem;

/**
 * Return the bounding rect for the UnboundMonomericAttachmentPointItem that
 * corresponds to the given parameters. Note that this graphics item doesn't
 * necessarily need to exist (and will *not* be created by this function).
 *
 * Parameter are the same as the UnboundMonomericAttachmentPointItem
 * constructor.
 */
SKETCHER_API QRectF get_bounding_rect_for_unbound_monomer_attachment_point_item(
    const UnboundAttachmentPoint& attachment_point,
    const AbstractMonomerItem* const parent_monomer, const Fonts& fonts);

/**
 * A graphics item for representing an unbound (available) attachment point
 * on a monomer. Draws a short line extending from the monomer with a text
 * label indicating the attachment point name.
 */
class SKETCHER_API UnboundMonomericAttachmentPointItem : public QGraphicsItem
{
  public:
    /**
     * @param attachment_point The attachment point to be represented by this
     * graphics item
     * @param parent_monomer The parent monomer graphics item
     * @param fonts The fonts to use for label rendering
     */
    UnboundMonomericAttachmentPointItem(
        const UnboundAttachmentPoint& attachment_point,
        AbstractMonomerItem* parent_monomer, const Fonts& fonts);

    enum {
        Type = static_cast<int>(
            GraphicsItemType::UNBOUND_MONOMERIC_ATTACHMENT_POINT_ITEM)
    };
    int type() const override;

    // QGraphicsItem overrides
    QRectF boundingRect() const override;
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

    /**
     * Set whether this attachment point indicator is active (black) or
     * inactive (gray).
     * @param active true for black coloring, false for gray
     */
    void setActive(bool active);

    /**
     * @return whether the given scene coordinates are within the attachment
     * point, but not within the parent monomer
     */
    bool withinHoverArea(const QPointF& scene_pos) const;

    /**
     * @return the attachment point represented by this graphics item
     */
    const UnboundAttachmentPoint& getAttachmentPoint() const;

  private:
    UnboundAttachmentPoint m_attachment_point;
    // this graphics item will never outlive the scene tool that created it, so
    // we store a pointer to the scene tool's fonts instead of making a copy
    const Fonts* m_fonts;
    QPointF m_line_end;
    QRectF m_label_rect;
    QString m_label_text;
    QRectF m_bounding_rect;
    QPainterPath m_hover_area;
    QPen m_line_pen;
    QBrush m_circle_brush{Qt::SolidPattern};
    bool m_is_active = false;

    /**
     * Update pen and brush colors based on active state.
     */
    void updateColors();
};

} // namespace sketcher
} // namespace schrodinger
