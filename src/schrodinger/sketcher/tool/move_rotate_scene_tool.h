#pragma once
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"
#include "schrodinger/sketcher/molviewer/rotation_item.h"
#include <GraphMol/Atom.h>
#include <Geometry/point.h>

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

  protected:
    void rotateRotationItem(QGraphicsSceneMouseEvent* const event);
    void rotate(QGraphicsSceneMouseEvent* const event, QPointF pivot_point,
                const std::vector<const RDKit::Atom*>& atoms_to_move = {});
    void translate(QGraphicsSceneMouseEvent* const event,
                   const std::vector<const RDKit::Atom*>& atoms_to_move = {});

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
     * Sets the atoms to be moved by this tool on a drag action. This is set
     * on mouse down by onDragStart and reset on mouse release by
     * onDragRelease and is used when a selection is present to distinguish
     * between actions that should be performed on the selection and those
     * that should be performed on the whole molecule
     */
    void setAtomsToMove(const std::vector<const RDKit::Atom*>& atoms_to_move);

    /** this variable is used to determine if the mouse is being performed in a
     * mouse drag (rotate, translate or none). It gets set in onDragStart and
     *reset to Action::NONE on mouse release
     */
    Action m_action = Action::NONE;

    RotationItem m_rotation_item = RotationItem();
    MoveSelectionItem m_move_selection_item = MoveSelectionItem();

    /**
     * compute the pivot point for a rotation. This is the selected atom that
     * shares a bond with the non-selected part of a molecule in case of the
     * selection of a terminal group. Otherwise it's the centroid of the
     * selection, or the centroid of the whole molecule if there's no selection
     */
    RDGeom::Point3D findPivotPointForRotation();

    std::vector<const RDKit::Atom*> m_atoms_to_move;
};

} // namespace sketcher
} // namespace schrodinger
