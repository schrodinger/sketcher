#pragma once

#include <QPen>
#include <QRectF>
#include <QString>

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

class SKETCHER_API AminoAcidItem : public AbstractAtomOrMonomerItem
{
  public:
    AminoAcidItem(const RDKit::Atom* atom, const Fonts& fonts,
                  QGraphicsItem* parent = nullptr);

    enum { Type = static_cast<int>(ItemType::AMINO_ACID) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

  protected:
    const Fonts& m_fonts;
    QRectF m_border_rect;
    QPen m_border_pen;
    QBrush m_border_brush = QBrush(Qt::BrushStyle::SolidPattern);
    QString m_main_label_text;
    QRectF m_main_label_rect;
    QFont m_main_label_font;
};

} // namespace sketcher
} // namespace schrodinger
