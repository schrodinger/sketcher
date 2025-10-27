#pragma once

#include <unordered_set>
#include <vector>

#include <QPen>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/Geometry/point.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/molviewer/rotation_item.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

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
class SKETCHER_API MoveRotateSceneTool : public StandardSceneToolBase
{
  public:
    MoveRotateSceneTool(Scene* scene, MolModel* mol_model);
    virtual void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onLeftButtonDragStart(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    updateColorsAfterBackgroundColorChange(bool has_dark_color_scheme) override;
    virtual void onSelectionChanged() override;
    virtual void onStructureUpdated() override;
    QPixmap createDefaultCursorPixmap() const override;

  protected:
    void updateGraphicsItemsDuringRotation(
        QGraphicsSceneMouseEvent* const event) override;
    void updateGraphicsItemsDuringTranslation(const QPointF& distance) override;

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
     * update the predictive highlighting when the mouse is moved to a new
     * position
     */
    void updatePredictiveHighlightingForMouseAt(const QPointF& pos);

    /** this variable is used to determine if the mouse is being performed in a
     * mouse drag (rotate, translate or none). It gets set in onDragStart and
     *reset to Action::NONE on mouse release
     */
    Action m_action = Action::NONE;

    RotationItem m_rotation_item = RotationItem();
    MoveSelectionItem m_move_selection_item = MoveSelectionItem();
};

} // namespace sketcher
} // namespace schrodinger
