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
 * The scene tool for right click mouse translation and context menus
 */
class SKETCHER_API TranslateSceneTool : public AbstractSceneTool
{
  public:
    TranslateSceneTool(Scene* scene, MolModel* mol_model);
    void onDragMove(QGraphicsSceneMouseEvent* const event) override;
    void onMouseClick(QGraphicsSceneMouseEvent* const event) override;
};

} // namespace sketcher
} // namespace schrodinger
