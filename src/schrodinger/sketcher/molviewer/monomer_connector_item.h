#pragma once

#include <QGraphicsItem>
#include <QLineF>
#include <QPen>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_bond_or_connector_item.h"

#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"

class QPainter;

namespace RDKit
{
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class AbstractMonomerItem;

/**
 * A Qt graphics item for representing a connector between two monomers in a
 * molviewer Scene.
 */
class SKETCHER_API MonomerConnectorItem : public AbstractBondOrConnectorItem
{
  public:
    /**
     * Note that this class does not take ownership of bond, so it must not be
     * destroyed while the graphics item is in use.
     *
     * @param bond The RDKit bond that this item should represent.
     *
     * @param start_monomer_item The graphics item for the starting monomer of
     * the bond
     *
     * @param end_monomer_item The graphics item for the ending monomer of the
     * bond
     *
     * @param is_secondary_connection Whether this graphics item represents the
     * secondary connection of the given bond. Secondary connections occur when
     * there is more than one connection between two monomers (e.g. neighboring
     * cysteines additionally joined by a disulfide bond). RDKit does not allow
     * more than one bond between two atoms, so a single bond object must
     * represent both connections.
     *
     * @param parent The Qt parent for this item.  See the QGraphicsItem
     * documentation for additional information.
     *
     * @pre bond != nullptr
     * @pre bond->hasOwningMol()
     */
    MonomerConnectorItem(const RDKit::Bond* bond,
                         const AbstractMonomerItem& start_monomer_item,
                         const AbstractMonomerItem& end_monomer_item,
                         const bool is_secondary_connection = false,
                         QGraphicsItem* parent = nullptr);

    enum { Type = static_cast<int>(ItemType::MONOMER_CONNECTOR) };
    int type() const override;

    // Overridden AbstractGraphicsItem method
    void updateCachedData() override;

    // Overridden QGraphicsItem methods
    void paint(QPainter* painter, const QStyleOptionGraphicsItem* option,
               QWidget* widget = nullptr) override;

    /**
     * @return whether this graphics item represents the secondary connection of
     * the bond
     */
    bool isSecondaryConnection() const;

  protected:
    QPen m_connector_pen;
    QPen m_arrowhead_pen;
    QBrush m_arrowhead_brush;
    QLineF m_connector_line;
    QPainterPath m_arrowhead_path;
    const AbstractMonomerItem& m_start_item;
    const AbstractMonomerItem& m_end_item;
    bool m_is_secondary_connection;
};

} // namespace sketcher
} // namespace schrodinger
