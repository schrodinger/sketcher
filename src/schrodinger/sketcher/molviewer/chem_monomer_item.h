#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API ChemMonomerItem : public AbstractAtomOrMonomerItem
{
  public:
    ChemMonomerItem(const RDKit::Atom* monomer, const Fonts& fonts,
                    QGraphicsItem* parent = nullptr);

    enum { Type = static_cast<int>(ItemType::CHEM_MONOMER) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

  protected:
    const Fonts& m_fonts;
};

} // namespace sketcher
} // namespace schrodinger
