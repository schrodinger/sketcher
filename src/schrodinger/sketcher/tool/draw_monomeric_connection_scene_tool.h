#pragma once

#include <optional>
#include <string>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/abstract_draw_monomer_or_monomeric_connection_scene_tool.h"

class QGraphicsItem;

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{

namespace rdkit_extensions
{
enum class ChainType;
} // namespace rdkit_extensions

namespace sketcher
{

/**
 * A scene tool that draws connections between monomers or monomeric attachment
 * points. Note that these connections use the standard attachments points, and
 * are therefore intended to represent a covalent or disulfide connection
 * between the monomers. The only exception to this is if a bond involves the
 * "pair" attachment point of a nucleic acids base, in which case it will
 * represent a hydrogen bond. To draw hydrogen bonds between other monomers, use
 * the DrawMonomericHBondSceneTool.
 */
class SKETCHER_API DrawMonomericConnectionSceneTool
    : public AbstractDrawMonomerOrMonomericConnectionSceneTool
{
  public:
    DrawMonomericConnectionSceneTool(const Fonts& fonts, Scene* scene,
                                     MolModel* mol_model);

    // Reimplemented AbstractSceneTool method
    void updateColorsAfterBackgroundColorChange(bool is_dark_mode) override;

  protected:
    QGraphicsItem* m_invalid_drag_item = nullptr;
    QPen m_invalid_drag_pen = QPen(CONNECTION_TOOL_INVALID_DRAG_COLOR,
                                   CONNECTION_TOOL_INVALID_DRAG_LINE_WIDTH);
    QColor m_unbound_ap_label_hover_color =
        CONNECTION_TOOL_AP_LABEL_HOVER_COLOR;

    // Reimplemented AbstractSceneTool method
    QPixmap createDefaultCursorPixmap() const override;

    // Reimplemented AbstractDrawMonomerOrMonomericConnectionSceneTool methods

    /**
     * Clicking never mutates a monomer for this tool.
     */
    bool clickShouldMutate(const RDKit::Atom* monomer,
                           const MonomerType monomer_type) const override;

    std::optional<HintFragmentMonomerInfo>
    getHintFragmentMonomerInfoForDragStart() override;

    std::optional<HintFragmentMonomerInfo> getHintFragmentMonomerInfoForDragEnd(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const DragEndInfo& drag_end_info) override;

    void updateHoveredUnboundAP(
        UnboundMonomericAttachmentPointItem* hovered_ap_item) override;

    DragEndInfo
    getDragEndInfoForNonMonomerPos(const QPointF& scene_pos) const override;

    void clearHintFragmentItem() override;

    void
    createHintFragmentItem(HintFragmentMonomerInfo& monomer_one_info,
                           HintFragmentMonomerInfo& monomer_two_info) override;

    void addDragStructureToMolModel(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const HintFragmentMonomerInfo& hint_end_monomer_info) override;

    UnboundMonomericAttachmentPointItem*
    getDefaultUnboundAttachmentPointForHoveredMonomer(
        const bool no_default_if_click_should_mutate) const override;
};

} // namespace sketcher
} // namespace schrodinger
