#pragma once

#include <optional>
#include <string>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/abstract_draw_monomeric_connection_scene_tool.h"

namespace schrodinger
{

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
    : public AbstractDrawMonomericConnectionSceneTool
{
  public:
    DrawMonomericConnectionSceneTool(const Fonts& fonts, Scene* scene,
                                     MolModel* mol_model);

  protected:
    // Reimplemented AbstractSceneTool method
    QPixmap createDefaultCursorPixmap() const override;

    std::optional<HintFragmentMonomerInfo> getHintFragmentMonomerInfoForDragEnd(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const DragEndInfo& drag_end_info) override;

    void updateHoveredUnboundAP(
        UnboundMonomericAttachmentPointItem* hovered_ap_item) override;

    UnboundMonomericAttachmentPointItem*
    getDefaultUnboundAttachmentPointForHoveredMonomer(
        const bool no_default_if_click_should_mutate) const override;
};

} // namespace sketcher
} // namespace schrodinger
