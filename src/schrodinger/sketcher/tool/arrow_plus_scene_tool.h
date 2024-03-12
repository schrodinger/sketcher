#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

namespace schrodinger
{
namespace sketcher
{

class Scene;
class MolModel;
enum class NonMolecularType;

/**
 * A scene tool for drawing arrows and plus signs
 */
class SKETCHER_API ArrowPlusSceneTool : public StandardSceneToolBase
{
  public:
    ArrowPlusSceneTool(const NonMolecularType& type, Scene* scene,
                       MolModel* mol_model);

    // reimplemented AbstractSceneTool method
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    QPixmap createDefaultCursorPixmap() const override;

  protected:
    NonMolecularType m_type;
};

} // namespace sketcher
} // namespace schrodinger
