#pragma once

#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A scene tool for adding/removing explicit Hs
 */
class SKETCHER_API ExplicitHsSceneTool : public StandardSceneToolBase
{
  public:
    ExplicitHsSceneTool(Scene* scene, MolModel* mol_model);

  protected:
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    QPixmap createDefaultCursorPixmap() const override;
};

} // namespace sketcher
} // namespace schrodinger
