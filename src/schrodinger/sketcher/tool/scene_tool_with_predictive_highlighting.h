#pragma once

#include <vector>

#include <QPointF>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{

class Scene;
class MolModel;

/**
 * A base class for scene tools that use predictive highlighting, where the atom
 * or bond under the mouse cursor is highlighted
 */
class SKETCHER_API SceneToolWithPredictiveHighlighting
    : public AbstractSceneTool
{
  public:
    SceneToolWithPredictiveHighlighting(Scene* scene, MolModel* mol_model);
    virtual ~SceneToolWithPredictiveHighlighting() = default;
    virtual void onMouseMove(QGraphicsSceneMouseEvent* event) override;
    virtual std::vector<QGraphicsItem*> getGraphicsItems() override;

  protected:
    PredictiveHighlightingItem m_predictive_highlighting_item =
        PredictiveHighlightingItem();
};

} // namespace sketcher
} // namespace schrodinger
