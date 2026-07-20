#pragma once

#include <optional>
#include <string>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/abstract_draw_monomeric_connection_scene_tool.h"

namespace schrodinger
{

namespace sketcher
{

/**
 * A scene tool that draws hydrogen bonds between monomers
 */
class SKETCHER_API DrawMonomericHBondSceneTool
    : public AbstractDrawMonomericConnectionSceneTool
{
  public:
    DrawMonomericHBondSceneTool(const Fonts& fonts, Scene* scene,
                                MolModel* mol_model);

  protected:
    // Reimplemented AbstractSceneTool method
    QPixmap createDefaultCursorPixmap() const override;

    void labelUnboundAttachmentPointsOnMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item,
        QGraphicsItemGroup& attachment_point_labels_group,
        std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items)
        override;

    // Reimplemented AbstractDrawMonomerOrMonomericConnectionSceneTool methods

    std::optional<HintFragmentMonomerInfo> getHintFragmentMonomerInfoForDragEnd(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const DragEndInfo& drag_end_info) override;

    void updateHoveredUnboundAP(
        UnboundMonomericAttachmentPointItem* hovered_ap_item) override;

    std::pair<DragEndInfo, AbstractMonomerItem*>
    getDragEndInfo(const QPointF& scene_pos) override;

    std::string
    getAPModelNameAtPosOverMonomer(const QPointF& scene_pos) const override;

    AbstractMonomerItem*
    getMonomerNear(const QPointF& scene_pos) const override;
};

} // namespace sketcher
} // namespace schrodinger
