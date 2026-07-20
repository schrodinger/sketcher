#pragma once

#include <optional>

#include <QColor>
#include <QPen>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/abstract_draw_monomer_or_monomeric_connection_scene_tool.h"

class QGraphicsItem;

namespace schrodinger
{
namespace sketcher
{

/**
 * A base class for scene tools that draw connections between existing monomers
 * without adding new monomers
 */
class SKETCHER_API AbstractDrawMonomericConnectionSceneTool
    : public AbstractDrawMonomerOrMonomericConnectionSceneTool
{
  protected:
    AbstractDrawMonomericConnectionSceneTool(const Fonts& fonts, Scene* scene,
                                             MolModel* mol_model);

    QGraphicsItem* m_invalid_drag_item = nullptr;
    QPen m_invalid_drag_pen = QPen(CONNECTION_TOOL_INVALID_DRAG_COLOR,
                                   CONNECTION_TOOL_INVALID_DRAG_LINE_WIDTH);
    QColor m_unbound_ap_label_hover_color =
        CONNECTION_TOOL_AP_LABEL_HOVER_COLOR;

    // Reimplemented AbstractSceneTool method
    void updateColorsAfterBackgroundColorChange(bool is_dark_mode) override;

    // Reimplemented AbstractDrawMonomerOrMonomericConnectionSceneTool methods
    bool clickShouldMutate(const RDKit::Atom* monomer,
                           MonomerType monomer_type) const override;

    std::optional<HintFragmentMonomerInfo>
    getHintFragmentMonomerInfoForDragStart() override;

    DragEndInfo
    getDragEndInfoForNonMonomerPos(const QPointF& scene_pos) const override;

    void clearHintFragmentItem() override;

    void
    createHintFragmentItem(HintFragmentMonomerInfo& monomer_one_info,
                           HintFragmentMonomerInfo& monomer_two_info) override;

    void addDragStructureToMolModel(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const HintFragmentMonomerInfo& hint_end_monomer_info) override;
};

} // namespace sketcher
} // namespace schrodinger
