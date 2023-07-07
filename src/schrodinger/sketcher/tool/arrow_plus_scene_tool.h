#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

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
class SKETCHER_API ArrowPlusSceneTool : public AbstractSceneTool
{
  public:
    ArrowPlusSceneTool(const NonMolecularType& type, Scene* scene,
                       MolModel* mol_model);

    // reimplemented AbstractSceneTool method
    void onMouseClick(QGraphicsSceneMouseEvent* const event) override;

  protected:
    NonMolecularType m_type;
};

} // namespace sketcher
} // namespace schrodinger
