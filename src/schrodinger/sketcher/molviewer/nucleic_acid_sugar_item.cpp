#include "schrodinger/sketcher/molviewer/nucleic_acid_sugar_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

NucleicAcidSugarItem::NucleicAcidSugarItem(const RDKit::Atom* monomer,
                                           const Fonts& fonts,
                                           QGraphicsItem* parent) :
    AbstractAtomOrMonomerItem(monomer, parent),
    m_fonts(fonts)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    updateCachedData();
}

int NucleicAcidSugarItem::type() const
{
    return Type;
}

void NucleicAcidSugarItem::updateCachedData()
{
    prepareGeometryChange();
}

void NucleicAcidSugarItem::paint(QPainter* painter,
                                 const QStyleOptionGraphicsItem* option,
                                 QWidget* widget)
{
}

} // namespace sketcher
} // namespace schrodinger
