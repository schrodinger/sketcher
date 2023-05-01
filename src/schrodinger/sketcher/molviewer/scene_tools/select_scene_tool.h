#pragma once

#include <memory>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/scene_tools/scene_tool_with_predictive_highlighting.h"
#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_items.h"

namespace schrodinger
{
namespace sketcher
{

class Scene;
class MolModel;
enum class SelectionTool;
enum class SelectMode;

/**
 * Instantiate and return the appropriate scene tool for the specified
 * SketcherModel tool.
 * @param selection_type The SketcherModel tool
 * @param scene The scene that this tool is for
 * @param mol_model The molecule displayed in scene
 * @return The newly instantiated scene tool
 */
std::shared_ptr<AbstractSceneTool>
get_select_scene_tool(SelectionTool selection_type, Scene* scene,
                      MolModel* mol_model);

/**
 * A base class for the lasso and marquee selection scene tools
 * @tparam T The graphics items class used to represent the user's selection as
 * it's being drawn
 */
template <typename T> class SelectSceneTool
    : public SceneToolWithPredictiveHighlighting
{
  public:
    SelectSceneTool(Scene* scene, MolModel* mol_model);
    virtual ~SelectSceneTool() = default;

    virtual void onDragStart(QGraphicsSceneMouseEvent* event) override;
    virtual void onDragMove(QGraphicsSceneMouseEvent* event) override;
    virtual void onDragRelease(QGraphicsSceneMouseEvent* event) override;

    virtual void onMouseClick(QGraphicsSceneMouseEvent* event) override;

    virtual std::vector<QGraphicsItem*> getGraphicsItems() override;

    /**
     * @return the MolModel select mode that corresponds to the user's current
     * keyboard modifiers, i.e., is Ctrl or Shift being held down.
     */
    SelectMode getSelectMode(QGraphicsSceneMouseEvent* event);

  protected:
    T m_select_item = T();
};

/**
 * The scene tool for lasso selection
 */
class LassoSelectSceneTool : public SelectSceneTool<LassoSelectionItem>
{
  public:
    LassoSelectSceneTool(Scene* scene, MolModel* mol_model);
    void onMousePress(QGraphicsSceneMouseEvent* event) override;
    void onMouseMove(QGraphicsSceneMouseEvent* event) override;
    void onDragRelease(QGraphicsSceneMouseEvent* event) override;
};

/**
 * The base class for marquee selection using a specified shape (rectangular or
 * elliptical)
 * @tparam T The graphics items class used to represent the user's selection as
 * it's being drawn
 */
template <typename T> class ShapeSelectSceneTool : public SelectSceneTool<T>
{
  public:
    ShapeSelectSceneTool(Scene* scene, MolModel* mol_model);
    void onDragMove(QGraphicsSceneMouseEvent* event) override;
};

typedef ShapeSelectSceneTool<RectSelectionItem> RectSelectSceneTool;
typedef ShapeSelectSceneTool<EllipseSelectionItem> EllipseSelectSceneTool;

} // namespace sketcher
} // namespace schrodinger
