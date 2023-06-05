#pragma once
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{
class Scene;
class MolModel;

/**
 * The scene tool for mid mouse rotation
 */
class SKETCHER_API RotateSceneTool : public AbstractSceneTool
{
  public:
    RotateSceneTool(Scene* scene, MolModel* mol_model);
    virtual void onDragMove(QGraphicsSceneMouseEvent* const event) override;
};

} // namespace sketcher
} // namespace schrodinger
