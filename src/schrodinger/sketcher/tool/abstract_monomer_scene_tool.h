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
 * Release the tool's references to attachment-point label items after the
 * parent monomer has just been destroyed by a structure update. Bound AP
 * labels in `group` were reparented to the group by addToGroup, so they
 * survive the monomer's destruction and must still be deleted normally.
 * Items in `unbound_ap_items` were Qt children of the now-gone monomer and
 * have already been destroyed by Qt — drop the dangling pointers without
 * deleting. Compare to clear_graphics_item_group_and_list (defined later),
 * which deletes both and is only safe while the parent monomer is alive.
 */
template <typename T> SKETCHER_API void
release_after_parent_destroyed(QGraphicsItemGroup& group,
                               std::vector<T*>& unbound_ap_items);

/**
 * Abstract base class for scene tools that work with monomers.
 * Provides common functionality for finding monomer items, managing hint
 * fragments, and drawing attachment point labels.
 */
class SKETCHER_API AbstractMonomerSceneTool : public StandardSceneToolBase
{
  public:
    AbstractMonomerSceneTool(const Fonts& fonts, Scene* scene,
                             MolModel* mol_model);
    virtual ~AbstractMonomerSceneTool();

    // Reimplemented StandardSceneToolBase methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    void updateColorsAfterBackgroundColorChange(bool is_dark_mode) override;
    virtual void onStructureUpdated() override;

  protected:
    Fonts m_fonts;
    MonomerHintFragmentItem* m_hint_fragment_item = nullptr;
    QGraphicsItemGroup m_attachment_point_labels_group;
    const UnboundMonomericAttachmentPointItem* m_hovered_ap_item = nullptr;
    std::vector<UnboundMonomericAttachmentPointItem*> m_unbound_ap_items;
    QColor m_unbound_ap_label_color = UNBOUND_AP_LABEL_COLOR;
    QColor m_bound_ap_label_color = BOUND_AP_LABEL_COLOR;
    const QGraphicsItem* m_hovered_item = nullptr;
    QColor m_monomer_background_color = LIGHT_BACKGROUND_COLOR;
    bool m_cursor_hint_shown = true;

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

    /**
     * Clear any existing attachment point labels and draw new ones for the
     * specified monomeric graphics item. If item is nullptr, then all labels
     * will be cleared and no new ones will be drawn.
     */
    void drawAttachmentPointLabelsFor(QGraphicsItem* const item);

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
     * Label all attachment points on the given hovered monomer
     */
    void labelAttachmentPointsOnHoveredMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item);

    /**
     * Label all attachment points on the given monomer, placing all newly
     * created graphics items into the specified group and list. (You probably
     * want to call either labelAttachmentPointsOnHoveredMonomer or
     * labelAttachmentPointsOnDragEndMonomer instead of this.)
     */
    void labelAttachmentPointsOnMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item,
        QGraphicsItemGroup& attachment_point_labels_group,
        std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items);

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
     * Show or hide the cursor hint
     */
    void setCursorHintShown(bool show);
};

} // namespace sketcher
} // namespace schrodinger
