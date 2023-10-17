#pragma once

#include <GraphMol/ROMol.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/tool/scene_tool_with_predictive_highlighting.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A scene tool for adding/removing explicit Hs
 */
class SKETCHER_API ExplicitHsSceneTool
    : public SceneToolWithPredictiveHighlighting
{
  public:
    ExplicitHsSceneTool(Scene* scene, MolModel* mol_model);

  protected:
    void onMouseClick(QGraphicsSceneMouseEvent* const event) override;
    QPixmap getCursorPixmap() const override;
};

} // namespace sketcher
} // namespace schrodinger
