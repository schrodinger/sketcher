#pragma once
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/move_rotate_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{
class Scene;
class MolModel;

/**
 * The scene tool for mouse rotation (middle button drag)
 */
class SKETCHER_API RotateSceneTool : public MoveRotateSceneTool
{
  public:
    RotateSceneTool(Scene* scene, MolModel* mol_model);
    void onDragStart(QGraphicsSceneMouseEvent* event) override;
    virtual void onDragMove(QGraphicsSceneMouseEvent* const event) override;
};

} // namespace sketcher
} // namespace schrodinger
