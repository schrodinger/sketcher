#include "schrodinger/sketcher/molviewer/abstract_bond_or_connector_item.h"

#include <QPointF>

namespace schrodinger
{
namespace sketcher
{

AbstractBondOrConnectorItem::AbstractBondOrConnectorItem(
    const RDKit::Bond* bond, QGraphicsItem* parent) :
    AbstractGraphicsItem(parent),
    m_bond(bond)
{
}

const RDKit::Bond* AbstractBondOrConnectorItem::getBond() const
{
    return m_bond;
}

QPointF AbstractBondOrConnectorItem::getMidpoint() const
{
    return m_midpoint;
}

} // namespace sketcher
} // namespace schrodinger
