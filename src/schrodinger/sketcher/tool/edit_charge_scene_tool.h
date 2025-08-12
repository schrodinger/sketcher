#pragma once

#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A scene tool for editing atom charges
 */
class SKETCHER_API EditChargeSceneTool : public StandardSceneToolBase
{
  public:
    EditChargeSceneTool(ChargeTool bond_tool, Scene* scene,
                        MolModel* mol_model);

  protected:
    // reimplemented AbstractSceneTool methods
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    QPixmap createDefaultCursorPixmap() const override;

    ChargeTool m_charge_tool;
};

} // namespace sketcher
} // namespace schrodinger
