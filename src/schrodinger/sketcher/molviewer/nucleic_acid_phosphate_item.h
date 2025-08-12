#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API NucleicAcidPhosphateItem : public AbstractMonomerItem
{
  public:
    NucleicAcidPhosphateItem(const RDKit::Atom* atom, const Fonts& fonts,
                             QGraphicsItem* parent = nullptr);

    enum { Type = static_cast<int>(ItemType::NUCLEIC_ACID_PHOSPHATE) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;
};

} // namespace sketcher
} // namespace schrodinger
