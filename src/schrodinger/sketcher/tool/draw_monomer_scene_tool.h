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

/**
 * A scene tools that draws a monomer
 */
class SKETCHER_API DrawMonomerSceneTool : public StandardSceneToolBase
{
  public:
    DrawMonomerSceneTool(const std::string& res_name,
                         const rdkit_extensions::ChainType chain_type,
                         const Fonts& fonts, Scene* scene, MolModel* mol_model);

    // Reimplemented AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;

  protected:
    std::string m_res_name;
    rdkit_extensions::ChainType m_chain_type;
    Fonts m_fonts;
    QGraphicsItemGroup m_attachment_point_labels_group;
    const QGraphicsItem* m_hovered_item = nullptr;

    QPixmap createDefaultCursorPixmap() const override;

    /**
     * Label all attachment points on the given monomer
     */
    void labelAttachmentPointsOnMonomer(const RDKit::Atom* const monomer);

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
     * Clear all attachment point labels drawn by this scene tool
     */
    void clearAttachmentPointsLabels();
};

} // namespace sketcher
} // namespace schrodinger
