#pragma once

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

class SKETCHER_API AminoAcidItem : public AbstractMonomerItem
{
  public:
    AminoAcidItem(const RDKit::Atom* atom, const Fonts& fonts,
                  QGraphicsItem* parent = nullptr);

    enum { Type = static_cast<int>(ItemType::AMINO_ACID) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;
};

} // namespace sketcher
} // namespace schrodinger
