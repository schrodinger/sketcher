#pragma once

#include <unordered_set>
#include <vector>

#include <QPen>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/Geometry/point.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/molviewer/rotation_item.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{
class Scene;
class MolModel;

class MoveSelectionItem : public QGraphicsRectItem
{
  public:
    MoveSelectionItem();
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
     * Draw circles centered at each of the given coordinates.  All previously
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

enum class Action { ROTATE, TRANSLATE, NONE };

/**
 * The scene tool for MOVE_ROTATE draw tool
 */
class SKETCHER_API MoveRotateSceneTool : public AbstractSceneTool
{
  public:
    MoveRotateSceneTool(Scene* scene, MolModel* mol_model);
    virtual void onDragStart(QGraphicsSceneMouseEvent* const event) override;
    virtual void onDragMove(QGraphicsSceneMouseEvent* const event) override;
    virtual void onDragRelease(QGraphicsSceneMouseEvent* const event) override;

    // TODO: remove this function, for testing only until context menus are
    // implemented
    virtual void onMouseClick(QGraphicsSceneMouseEvent* const event) override;

    virtual void onStructureUpdated() override;
    QPixmap getCursorPixmap() const override;

  protected:
    void rotateRotationItem(QGraphicsSceneMouseEvent* const event);
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

    virtual std::vector<QGraphicsItem*> getGraphicsItems() override;

    /**
     * make sure the rotation item is centered on the pivot point and is
     * consistent with the current molecule and selection. This function is
     * called by getGraphicsItem, right before the item is returned to be
     * displayed by the scene
     */
    void updateRotationItem();

    /**
     * make sure the move selection item depicts the outline of the current
     * selection. This is updated on every drag move, since both rotation and
     * translation can change the selection's coordinates
     */
    void updateMoveSelectionItem();

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

    /** this variable is used to determine if the mouse is being performed in a
     * mouse drag (rotate, translate or none). It gets set in onDragStart and
     *reset to Action::NONE on mouse release
     */
    Action m_action = Action::NONE;

    RotationItem m_rotation_item = RotationItem();
    MoveSelectionItem m_move_selection_item = MoveSelectionItem();
    MergeHintItem m_merge_hint_item = MergeHintItem();

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
 * Contains the implementation for MoveRotateSceneTool::getOverlappingAtomIdxs
 * for unit testing purposes.
 */
SKETCHER_API std::vector<std::pair<unsigned int, unsigned int>>
get_overlapping_atom_idxs(
    const RDKit::ROMol* const mol,
    const std::unordered_set<const RDKit::Atom*> atoms_to_move);

} // namespace sketcher
} // namespace schrodinger
