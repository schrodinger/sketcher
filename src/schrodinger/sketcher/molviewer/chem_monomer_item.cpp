#include "schrodinger/sketcher/molviewer/chem_monomer_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

ChemMonomerItem::ChemMonomerItem(const RDKit::Atom* monomer, const Fonts& fonts,
                                 QGraphicsItem* parent) :
    AbstractAtomOrMonomerItem(monomer, parent),
    m_fonts(fonts)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    updateCachedData();
}

int ChemMonomerItem::type() const
{
    return Type;
}

void ChemMonomerItem::updateCachedData()
{
    prepareGeometryChange();
}

void ChemMonomerItem::paint(QPainter* painter,
                            const QStyleOptionGraphicsItem* option,
                            QWidget* widget)
{
}

} // namespace sketcher
} // namespace schrodinger
