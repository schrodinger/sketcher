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

/**
 * A graphics item for representing a CHEM monomer (i.e. a monomer that isn't an
 * amino acid or a nucleic acid) as a rectangle
 */
class SKETCHER_API ChemMonomerItem : public AbstractMonomerItem
{
  public:
    ChemMonomerItem(const RDKit::Atom* monomer, const Fonts& fonts,
                    QGraphicsItem* parent = nullptr);

    enum { Type = static_cast<int>(ItemType::CHEM_MONOMER) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;
};

} // namespace sketcher
} // namespace schrodinger
