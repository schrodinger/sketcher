#include "schrodinger/sketcher/molviewer/scene_tools/select_scene_tool.h"

#include <QList>
#include <QGraphicsItem>

#include "schrodinger/sketcher/qt_utils.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/mol_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"

namespace schrodinger
{
namespace sketcher
{

std::shared_ptr<AbstractSceneTool>
get_select_scene_tool(SelectionTool selection_type, Scene* scene,
                      MolModel* mol_model)
{
    if (selection_type == SelectionTool::RECTANGLE) {
        return std::make_shared<RectSelectSceneTool>(scene, mol_model);
    } else if (selection_type == SelectionTool::ELLIPSE) {
        return std::make_shared<EllipseSelectSceneTool>(scene, mol_model);
    } else { // selection_type == SelectionTool::LASSO
        return std::make_shared<LassoSelectSceneTool>(scene, mol_model);
    }
}

template <typename T>
SelectSceneTool<T>::SelectSceneTool(Scene* scene, MolModel* mol_model) :
    SceneToolWithPredictiveHighlighting(scene, mol_model)
{
    m_select_item.setVisible(false);
}

template <typename T>
void SelectSceneTool<T>::onDragStart(QGraphicsSceneMouseEvent* event)
{
    SceneToolWithPredictiveHighlighting::onDragStart(event);
    m_select_item.setVisible(true);
}

template <typename T>
void SelectSceneTool<T>::onDragMove(QGraphicsSceneMouseEvent* event)
{
    SceneToolWithPredictiveHighlighting::onDragMove(event);
    QList<QGraphicsItem*> items = m_scene->collidingItems(&m_select_item);
    m_predictive_highlighting_item.highlightItems(items);
}

template <typename T>
void SelectSceneTool<T>::onDragRelease(QGraphicsSceneMouseEvent* event)
{
    SceneToolWithPredictiveHighlighting::onDragRelease(event);
    m_select_item.setVisible(false);
    m_predictive_highlighting_item.clearHighlightingPath();

    auto select_mode = getSelectMode(event);
    QList<QGraphicsItem*> items = m_scene->collidingItems(&m_select_item);
    m_scene->selectGraphicsItems(items, select_mode);
}

template <typename T>
SelectMode SelectSceneTool<T>::getSelectMode(QGraphicsSceneMouseEvent* event)
{
    auto modifiers = event->modifiers();
    if (modifiers & Qt::ControlModifier) {
        return SelectMode::TOGGLE;
    } else if (modifiers & Qt::ShiftModifier) {
        return SelectMode::SELECT;
    } else {
        return SelectMode::SELECT_ONLY;
    }
}

template <typename T>
void SelectSceneTool<T>::onMouseClick(QGraphicsSceneMouseEvent* event)
{
    SceneToolWithPredictiveHighlighting::onMouseClick(event);
    auto select_mode = getSelectMode(event);
    QGraphicsItem* item = m_scene->itemAt(event->scenePos(), QTransform());
    if (item != nullptr) {
        m_scene->selectGraphicsItems({item}, select_mode);
    } else {
        m_scene->selectGraphicsItems({}, select_mode);
    }
}

template <typename T>
std::vector<QGraphicsItem*> SelectSceneTool<T>::getGraphicsItems()
{
    auto items = SceneToolWithPredictiveHighlighting::getGraphicsItems();
    items.push_back(&m_select_item);
    return items;
}

LassoSelectSceneTool::LassoSelectSceneTool(Scene* scene, MolModel* mol_model) :
    SelectSceneTool(scene, mol_model)
{
}

void LassoSelectSceneTool::onMousePress(QGraphicsSceneMouseEvent* event)
{
    SelectSceneTool<LassoSelectionItem>::onMousePress(event);
    m_select_item.clearPath();
    m_select_item.addPoint(event->scenePos());
}

void LassoSelectSceneTool::onMouseMove(QGraphicsSceneMouseEvent* event)
{
    SelectSceneTool<LassoSelectionItem>::onMouseMove(event);
    // add a point to the path even if we haven't started a drag yet.  (The drag
    // doesn't start until the mouse has moved at least
    // QApplication::startDragDistance() pixels, but we want the path to include
    // all of the cursor movement before that once the path is shown.)
    if (m_mouse_pressed) {
        m_select_item.addPoint(event->scenePos());
    }
}

void LassoSelectSceneTool::onDragRelease(QGraphicsSceneMouseEvent* event)
{
    SelectSceneTool<LassoSelectionItem>::onDragRelease(event);
    m_select_item.clearPath();
}

template <typename T>
ShapeSelectSceneTool<T>::ShapeSelectSceneTool(Scene* scene,
                                              MolModel* mol_model) :
    SelectSceneTool<T>(scene, mol_model)
{
}

template <typename T>
void ShapeSelectSceneTool<T>::onDragMove(QGraphicsSceneMouseEvent* event)
{
    QRectF rect =
        rect_for_points(this->m_mouse_press_scene_pos, event->scenePos());
    this->m_select_item.setRect(rect);
    SelectSceneTool<T>::onDragMove(event);
}

} // namespace sketcher
} // namespace schrodinger
