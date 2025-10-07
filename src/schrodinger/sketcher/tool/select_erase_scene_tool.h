#pragma once

#include <memory>
#include <tuple>
#include <unordered_set>
#include <utility>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_items.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class Scene;
class MolModel;
class NonMolecularObject;
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
SKETCHER_API std::shared_ptr<AbstractSceneTool>
get_select_scene_tool(SelectionTool selection_type, Scene* scene,
                      MolModel* mol_model);

/**
 * A base class for the lasso and marquee selection scene tools, as well as the
 * erase scene tool (which inherits from RectSelectSceneTool)
 * @tparam T The graphics items class used to represent the user's selection as
 * it's being drawn
 */
template <typename T> class SKETCHER_API SelectSceneTool
    : public StandardSceneToolBase
{
  public:
    SelectSceneTool(Scene* scene, MolModel* mol_model);

    virtual void
    onLeftButtonDragStart(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onRightButtonClick(QGraphicsSceneMouseEvent* const event) override;

    virtual void
    onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    virtual void
    onLeftButtonDoubleClick(QGraphicsSceneMouseEvent* const event) override;

    virtual void
    updateColorsAfterBackgroundColorChange(bool is_dark_mode) override;

    virtual std::vector<QGraphicsItem*> getGraphicsItems() override;

  protected:
    T m_select_item = T();

    /**
     * @return the MolModel select mode that corresponds to the user's current
     * keyboard modifiers, i.e., is Ctrl or Shift being held down.
     */
    SelectMode getSelectMode(QGraphicsSceneMouseEvent* const event) const;

    /**
     * Carry out the MolModel selection after the user has selected item(s) in
     * the Scene.  Note that this method is overridden in the EraseSceneTool
     * subclass to erase rather than select.
     * @param items The graphics items that were selected in the Scene
     * @param event The event that triggered the selection
     */
    virtual void onSelectionMade(const QList<QGraphicsItem*>& items,
                                 QGraphicsSceneMouseEvent* event);
};

/**
 * The scene tool for lasso selection
 */
class SKETCHER_API LassoSelectSceneTool
    : public SelectSceneTool<LassoSelectionItem>
{
  public:
    LassoSelectSceneTool(Scene* scene, MolModel* mol_model);
    void onLeftButtonPress(QGraphicsSceneMouseEvent* const event) override;
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;
    QPixmap createDefaultCursorPixmap() const override;
};

/**
 * The base class for marquee selection using a specified shape (rectangular or
 * elliptical)
 * @tparam T The graphics items class used to represent the user's selection as
 * it's being drawn
 */
template <typename T> class SKETCHER_API ShapeSelectSceneTool
    : public SelectSceneTool<T>
{
  public:
    ShapeSelectSceneTool(Scene* scene, MolModel* mol_model);
    void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    QPixmap createDefaultCursorPixmap() const override;
};

typedef ShapeSelectSceneTool<RectSelectionItem> RectSelectSceneTool;
typedef ShapeSelectSceneTool<EllipseSelectionItem> EllipseSelectSceneTool;

// We define EraseSceneTool in the select scene tool file to avoid issues with
// inheriting from a templated class across compilation units
class SKETCHER_API EraseSceneTool : public RectSelectSceneTool
{
  public:
    EraseSceneTool(Scene* scene, MolModel* mol_model);
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDoubleClick(QGraphicsSceneMouseEvent* const event) override;
    void onRightButtonClick(QGraphicsSceneMouseEvent* const event) override;
    QPixmap createDefaultCursorPixmap() const override;

  protected:
    void onSelectionMade(const QList<QGraphicsItem*>& items,
                         QGraphicsSceneMouseEvent* const event) override;
};

} // namespace sketcher
} // namespace schrodinger
