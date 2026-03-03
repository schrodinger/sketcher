#pragma once

#include <string>
#include <tuple>
#include <utility>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

namespace schrodinger
{

namespace rdkit_extensions
{
enum class ChainType;
} // namespace rdkit_extensions

namespace sketcher
{

enum class MonomerType;
class UnboundMonomericAttachmentPointItem;

/**
 * @return the attachment point name to a QString after converting apostrophes
 * to Unicode primes.
 */
SKETCHER_API QString prep_attachment_point_name(const std::string& name);

/**
 * Position the given rectangle to label a monomer's attachment point
 * @param ap_label_rect The rectangle to position. It should already be sized
 * correctly for the attachment point label.
 * @param monomer_coords The coordinates of the monomer being labeled
 * @param bound_coords The coordinates of the other monomer involved in the bond
 */
SKETCHER_API void position_ap_label_rect(QRectF& ap_label_rect,
                                         const QPointF& monomer_coords,
                                         const QPointF& bound_coords);

/**
 * Return the default unbound attachment point; that is, the attachment point
 * that should be selected when the user hovers over a monomer.
 * @param hovered_type The type of monomer being hovered over
 * @param tool_type  The type of monomer that would be drawn by the active scene
 * tool
 * @param unbound_ap_items A list of all graphics items representing unbound
 * attachment points of the hovered monomer
 */
SKETCHER_API UnboundMonomericAttachmentPointItem* get_default_attachment_point(
    const MonomerType hovered_type, const MonomerType tool_type,
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items);

/**
 * A scene tools that draws a monomer
 */
class SKETCHER_API DrawMonomerSceneTool : public StandardSceneToolBase
{
  public:
    DrawMonomerSceneTool(const std::string& res_name,
                         const rdkit_extensions::ChainType chain_type,
                         const Fonts& fonts, Scene* scene, MolModel* mol_model);
    virtual ~DrawMonomerSceneTool();

    // Reimplemented AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;

  protected:
    std::string m_res_name;
    rdkit_extensions::ChainType m_chain_type;
    Fonts m_fonts;
    MonomerType m_monomer_type;
    QGraphicsItemGroup m_attachment_point_labels_group;
    const QGraphicsItem* m_hovered_item = nullptr;
    std::vector<UnboundMonomericAttachmentPointItem*> m_unbound_ap_items;

    QPixmap createDefaultCursorPixmap() const override;

    /**
     * Label all attachment points on the given monomer
     */
    void
    labelAttachmentPointsOnMonomer(const RDKit::Atom* const monomer,
                                   AbstractMonomerItem* const monomer_item);

    /**
     * Label both attachment points for the given monomeric connector
     * @param connector The monomeric connector to label
     * @param is_secondary_connection Whether this method should label the
     * secondary connection of the bond instead of the primary connection.
     * Secondary connections occur when a single RDKit::Bond* represents
     * multiple connections, e.g. two neighboring cysteines that are disulfide
     * bonded to each other.
     */
    void labelAttachmentPointsOnConnector(const RDKit::Bond* const connector,
                                          const bool is_secondary_connection);

    /**
     * Label the specified attachment point
     * @param monomer The monomer containing the attachment point to label
     * @param bound_monomer The other monomer involved in the connection
     * @param is_secondary_connection Whether we should label the secondary
     * connection of the bond instead of the primary connection.
     * @param label The text to display in the attachment point label
     */
    void labelBoundAttachmentPoint(const RDKit::Atom* const monomer,
                                   const RDKit::Atom* const bound_monomer,
                                   const bool is_secondary_connection,
                                   const std::string& label);

    /**
     * Place an attachment point label at the center of the given monomeric
     * connector.  This is used to label nucleic acid base pairs, as both
     * attachment points have the same label ("pair") and the bond is typically
     * too small to fit two labels.
     */
    void labelCenterOfConnector(const RDKit::Atom* const begin_monomer,
                                const RDKit::Atom* const end_monomer,
                                const QString& label);

    void addAttachmentPointLabel(const QString& label,
                                 const QRectF& label_rect);

    /**
     * Clear all attachment point labels drawn by this scene tool
     */
    void clearAttachmentPointsLabels();

    /**
     * Return the top graphics item representing a monomer or monomeric
     * connector at the given coordinates. Note that this method will consider
     * the cursor to be over a monomer if the cursor is over an unbound
     * attachment point belonging to that monomer, or if the cursor *would be*
     * over an unbound attachment point once it's drawn.
     * @param scene_pos The position in Scene coordinates
     */
    QGraphicsItem* getTopMonomericItemAt(const QPointF& scene_pos);

    /**
     * Clear any existing attachment point labels and draw new ones for the
     * specified monomeric graphics item. If item is nullptr, then all labels
     * will be cleared and no new ones will be drawn.
     */
    void drawAttachmentPointLabelsFor(QGraphicsItem* const item);

    /**
     * Return the unbound attachment point graphics item that should be "active"
     * (i.e. that we'd drawn a connection for if the user clicked) for the given
     * coordinates. If the coordinates are over an unbound attachment point
     * item, that item will be returned. If the coordinates are over the
     * monomer, then the default unbound attachment point will be returned,
     * assuming one exists for the current tool. If no default attachment point
     * exists (e.g. if a click would mutate the monomer or if the tool and the
     * hovered monomer are of different molecule types), then nullptr will be
     * returned.
     *
     * @note This method only returns accurate results for the currently hovered
     * monomer, and assumes that the unbound attachment point graphics items
     * have already been drawn for this monomer.
     */
    UnboundMonomericAttachmentPointItem*
    getUnboundAttachmentPointAt(const QPointF& scene_pos);
};

} // namespace sketcher
} // namespace schrodinger
