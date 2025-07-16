#include "schrodinger/sketcher/molviewer/nucleic_acid_phosphate_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

NucleicAcidPhosphateItem::NucleicAcidPhosphateItem(const RDKit::Atom* monomer,
                                                   const Fonts& fonts,
                                                   QGraphicsItem* parent) :
    AbstractAtomOrMonomerItem(monomer, parent),
    m_fonts(fonts)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    updateCachedData();
}

int NucleicAcidPhosphateItem::type() const
{
    return Type;
}

void NucleicAcidPhosphateItem::updateCachedData()
{
    prepareGeometryChange();
}

void NucleicAcidPhosphateItem::paint(QPainter* painter,
                                     const QStyleOptionGraphicsItem* option,
                                     QWidget* widget)
{
}

} // namespace sketcher
} // namespace schrodinger
