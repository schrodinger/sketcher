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

class SKETCHER_API NucleicAcidSugarItem : public AbstractAtomOrMonomerItem
{
  public:
    NucleicAcidSugarItem(const RDKit::Atom* atom, const Fonts& fonts,
                         QGraphicsItem* parent = nullptr);

    enum { Type = static_cast<int>(ItemType::NUCLEIC_ACID_SUGAR) };
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
