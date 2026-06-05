#pragma once

#include <QColor>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
}

namespace schrodinger
{
namespace rdkit_extensions
{
enum class Direction;
}

namespace sketcher
{

class MonomerHintFragmentItem;
class MonomerConnectorItem;
class AbstractMonomerItem;
class UnboundMonomericAttachmentPointItem;

/**
 * @return coordinates that are roughly BOND_LENGTH units away from the given
 * monomer in the specified direction. (Monomers are laid out in a grid-like
 * pattern, so diagonal directions will lead to bonds that are slightly longer
 * than BOND_LENGTH.)
 */
SKETCHER_API RDGeom::Point3D
get_default_coords_for_bound_monomer(const RDGeom::Point3D monomer_pos,
                                     const rdkit_extensions::Direction dir);

/**
 * Abstract base class for scene tools that draw monomers. Provides common
 * functionality for drawing unbound attachment point labels and managing hint
 * fragments.
 */
class SKETCHER_API AbstractMonomerSceneTool : public StandardSceneToolBase
{
  public:
    AbstractMonomerSceneTool(const Fonts& fonts, Scene* scene,
                             MolModel* mol_model);
    virtual ~AbstractMonomerSceneTool();
    // Reimplemented StandardSceneToolBase methods
    void onStructureUpdated() override;
    void updateColorsAfterBackgroundColorChange(bool is_dark_mode) override;

  protected:
    MonomerHintFragmentItem* m_hint_fragment_item = nullptr;
    const UnboundMonomericAttachmentPointItem* m_hovered_ap_item = nullptr;
    std::vector<UnboundMonomericAttachmentPointItem*> m_unbound_ap_items;
    QColor m_unbound_ap_label_color = UNBOUND_AP_LABEL_COLOR;
    QColor m_monomer_background_color = LIGHT_BACKGROUND_COLOR;
    bool m_cursor_hint_shown = true;

    // we override this method so that hovering over an unbound attachment point
    // label counts as hovering over the monomer
    QGraphicsItem*
    getTopMonomericItemAt(const QPointF& scene_pos) const override;

    /**
     * Return the top graphics item representing a monomer at the given
     * coordinates. If there's no such graphics item, return nullptr.  Note that
     * this method will consider the cursor to be over a monomer if the cursor
     * is over an unbound attachment point belonging to that monomer, or if the
     * cursor *would be* over an unbound attachment point once it's drawn.
     * @param scene_pos The position in Scene coordinates
     */
    AbstractMonomerItem* getTopMonomerItemAt(const QPointF& scene_pos) const;

    /**
     * Return the top graphics item representing a monomeric connector at the
     * given coordinates. If there's no such graphics item, return nullptr.
     * @param scene_pos The position in Scene coordinates
     */
    MonomerConnectorItem*
    getTopMonomerConnectorItemAt(const QPointF& scene_pos) const;

    void
    drawMonomericAttachmentPointLabelsFor(QGraphicsItem* const item) override;

    /**
     * Remove the hint fragment item from the scene, then destroy it as well as
     * the RDKit structure object it represents.
     */
    void clearHintFragmentItem();

    /**
     * Determine whether the specified graphics item is part of this scene
     * tool's hint fragment
     */
    bool isItemPartOfHintFragment(QGraphicsItem* item) const;

    /**
     * Clear all attachment point labels drawn on the hovered monomer. Also
     * clear the hint fragment item and the associated structure.
     */
    void clearAttachmentPointsLabelsAndHintFragmentItem();
    /**
     * Label all unbound attachment points on the given monomer, placing all
     * newly created graphics items into the specified group and list.
     */
    void labelUnboundAttachmentPointsOnMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item,
        QGraphicsItemGroup& attachment_point_labels_group,
        std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items);

    /**
     * Show or hide the cursor hint
     */
    void setCursorHintShown(bool show);
};

void delete_all_unbound_ap_items(
    std::vector<UnboundMonomericAttachmentPointItem*>& items_list);

} // namespace sketcher
} // namespace schrodinger
