#pragma once

#include <optional>
#include <string>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/abstract_draw_monomer_or_monomeric_connection_scene_tool.h"

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

struct HintFragmentMonomerInfo;

/**
 * A scene tools that draws a monomer
 */
class SKETCHER_API DrawMonomerSceneTool
    : public AbstractDrawMonomerOrMonomericConnectionSceneTool
{
  public:
    DrawMonomerSceneTool(const std::string& res_name,
                         const rdkit_extensions::ChainType chain_type,
                         const Fonts& fonts, Scene* scene, MolModel* mol_model);

    // Reimplemented AbstractSceneTool method
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;

  protected:
    // Reimplemented AbstractSceneTool method
    QPixmap createDefaultCursorPixmap() const override;

    /**
     * Draw a hint structure showing a monomer bound to the specified attachment
     * point.
     */
    void
    drawBoundMonomerHintFor(UnboundMonomericAttachmentPointItem* const ap_item);

    // Reimplemented AbstractDrawMonomerOrMonomericConnectionSceneTool methods

    /**
     * Whether clicking on the specified monomer should mutate that monomer. We
     * only mutate monomers that are the same monomer type as this tool but have
     * a different residue name.
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
};

} // namespace sketcher
} // namespace schrodinger
