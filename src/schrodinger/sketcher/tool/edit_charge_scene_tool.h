#pragma once

#include <GraphMol/ROMol.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/tool/scene_tool_with_predictive_highlighting.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A scene tool for editing atom charges
 */
class SKETCHER_API EditChargeSceneTool
    : public SceneToolWithPredictiveHighlighting
{
  public:
    EditChargeSceneTool(ChargeTool bond_tool, Scene* scene,
                        MolModel* mol_model);

  protected:
    void onMouseClick(QGraphicsSceneMouseEvent* const event) override;

    ChargeTool m_charge_tool;
};

} // namespace sketcher
} // namespace schrodinger
