#pragma once

#include <unordered_set>
#include <vector>

#include <QPen>
#include <QPixmap>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/Geometry/point.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/rotation_item.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{
class Scene;
class MolModel;

/**
 * A graphics item that shows the rotation angle's value
 */
class AngleTextItem : public QGraphicsSimpleTextItem
{
  public:
    AngleTextItem();
    void centerOn(const QPointF& center);
};

/**
 * A graphics item group that displays circles around all pairs of atoms that
 * will be merged at the end of the drag.
 */
class MergeHintItem : public QGraphicsItemGroup
{
  public:
    MergeHintItem(QGraphicsItem* parent = nullptr);

    /**
     * Draw circles centered at each of the given coordinates. All previously
     * drawn circles will be cleared.
     */
    void setCoordinates(std::vector<QPointF> centers);

    /**
     * Clear all circles
     */
    void clear();

  protected:
    QPen m_circle_pen;
};

/**
 * A scene tool that contains the behavior shared between all (or at least most)
 * scene tools. With this scene tool, middle-button drags will rotate the
 * molecule, right-button drags will translate the molecule, and right-clicks
 * will open the context menu. Subclasses can additionally customize by setting
 * m_highlight_types (for predictive highlighting) or m_allow_context_menu (to
 * disable the context menu).
 */
class SKETCHER_API StandardSceneToolBase : public AbstractSceneTool
{
  public:
    StandardSceneToolBase(Scene* scene, MolModel* mol_model);
    virtual std::vector<QGraphicsItem*> getGraphicsItems() override;

  protected:
    /**
     * Whether this scene tool should open the context menu on a right-click.
     * Subclasses may change this to false to disable this behavior.
     */
    bool m_allow_context_menu = true;

    /**
     * The graphic item types that should be predictively highlighted. Default
     * behavior is NONE, i.e., no predictive highlighting. Subclasses may
     * change this to enable predictive highlighting.
     */
    InteractiveItemFlagType m_highlight_types = InteractiveItemFlag::NONE;

    PredictiveHighlightingItem m_predictive_highlighting_item =
        PredictiveHighlightingItem();
    QPixmap m_rotate_cursor_hint = cursor_hint_from_svg(":/icons/rotate.svg");
    QPixmap m_translate_cursor_hint =
        cursor_hint_from_svg(":/icons/translate.svg");

    virtual void
    onRightButtonClick(QGraphicsSceneMouseEvent* const event) override;
    virtual void onMouseMove(QGraphicsSceneMouseEvent* const event) override;

    virtual void
    onMiddleButtonDragStart(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onMiddleButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onMiddleButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;

    virtual void
    onRightButtonDragStart(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onRightButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onRightButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;

    virtual void
    updateColorsAfterBackgroundColorChange(bool is_dark_mode) override;

    void finishDrag();

    void
    rotate(QGraphicsSceneMouseEvent* const event, QPointF pivot_point,
           const std::unordered_set<const RDKit::Atom*>& atoms_to_move = {},
           const std::unordered_set<const NonMolecularObject*>&
               non_mol_objs_to_move = {});
    void
    translate(QGraphicsSceneMouseEvent* const event,
              const std::unordered_set<const RDKit::Atom*>& atoms_to_move = {},
              const std::unordered_set<const NonMolecularObject*>&
                  non_mol_objs_to_move = {});

    // These methods are no-ops in this class, but they're overridden in
    // MoveRotateSceneTool to update the rotation arm graphics item
    virtual void
    updateGraphicsItemsDuringRotation(QGraphicsSceneMouseEvent* const event);
    virtual void updateGraphicsItemsDuringTranslation(const QPointF& distance);

    /**
     * Update the merge hint item so it draws circles around all currently
     * overlapping atoms.
     */
    void updateMergeHintItem();

    /**
     * Sets the atoms to be moved by this tool on a drag action. This is set
     * on mouse down by onDragStart and reset on mouse release by
     * onDragRelease and is used when a selection is present to distinguish
     * between actions that should be performed on the selection and those
     * that should be performed on the whole molecule
     */
    void setObjectsToMove(
        const std::unordered_set<const RDKit::Atom*>& atoms_to_move,
        const std::unordered_set<const NonMolecularObject*>& non_mol_objects);

    /**
     * if a selection is present, set the selected objects to be moved by this
     * tool.
     */
    void setCurrentSelectionAsObjectsToMove();

    MergeHintItem m_merge_hint_item = MergeHintItem();
    AngleTextItem m_angle_text_item = AngleTextItem();

    /**
     * compute the pivot point for a rotation. This is the selected atom that
     * shares a bond with the non-selected part of a molecule in case of the
     * selection of a terminal group. Otherwise it's the centroid of the
     * selection, or the centroid of the whole molecule if there's no selection
     */
    RDGeom::Point3D findPivotPointForRotation();

    /**
     * @return the atom indices for all currently overlapping pairs of atoms.
     * The atom being moved will always be the first element of the pair.
     *
     * Two atoms are considered to be overlapping if
     *   - one atom is being moved and one atom is stationary
     *   - they are from different fragments
     *   - they are within MAX_DIST_FOR_DRAG_MERGE distance of each other.
     *
     * If an atom to move can potentially overlap with multiple stationary atoms
     * (e.g. there are two stationary atoms MAX_DIST_FOR_DRAG_MERGE from each
     * other), it will only be considered overlapping with the closest
     * stationary atom.
     */
    std::vector<std::pair<unsigned int, unsigned int>> getOverlappingAtomIdxs();

    /**
     * Merge all overlapping atoms at the end of a rotation or translation.
     */
    void mergeOverlappingAtoms();

    std::unordered_set<const RDKit::Atom*> m_atoms_to_move;
    std::unordered_set<const NonMolecularObject*> m_non_mol_objs_to_move;
};

/**
 * Contains the implementation for StandardSceneToolBase::getOverlappingAtomIdxs
 * for unit testing purposes.
 */
SKETCHER_API std::vector<std::pair<unsigned int, unsigned int>>
get_overlapping_atom_idxs(
    const RDKit::ROMol* const mol,
    const std::unordered_set<const RDKit::Atom*> atoms_to_move);

} // namespace sketcher
} // namespace schrodinger
