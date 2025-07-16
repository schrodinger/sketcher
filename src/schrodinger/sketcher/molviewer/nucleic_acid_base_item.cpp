#include "schrodinger/sketcher/molviewer/nucleic_acid_base_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

NucleicAcidBaseItem::NucleicAcidBaseItem(const RDKit::Atom* monomer,
                                         const Fonts& fonts,
                                         QGraphicsItem* parent) :
    AbstractAtomOrMonomerItem(monomer, parent),
    m_fonts(fonts)
{
    setZValue(static_cast<qreal>(ZOrder::MONOMER));
    updateCachedData();
}

int NucleicAcidBaseItem::type() const
{
    return Type;
}

void NucleicAcidBaseItem::updateCachedData()
{
    prepareGeometryChange();
}

void NucleicAcidBaseItem::paint(QPainter* painter,
                                const QStyleOptionGraphicsItem* option,
                                QWidget* widget)
{
}

} // namespace sketcher
} // namespace schrodinger
