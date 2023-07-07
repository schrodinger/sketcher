#pragma once

#include <QPainterPath>
#include <QPen>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API NonMolecularItem : public AbstractGraphicsItem
{
  public:
    NonMolecularItem(const NonMolecularObject* const non_molecular_object,
                     const QColor& color);

    // Type and type() are required for qgraphicsitem_cast.  Note that this enum
    // implementation is based on the sample code from the QGraphicsItem::Type
    // documentation.
    enum { Type = static_cast<int>(ItemType::NON_MOLECULAR) };
    int type() const override;

    /**
     * @return the non-molecular object associated with this item
     */
    const NonMolecularObject* getNonMolecularObject() const;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

  protected:
    const NonMolecularObject* m_non_molecular_object;
    QPainterPath m_paint_path;
    QPen m_pen;
};

} // namespace sketcher
} // namespace schrodinger
