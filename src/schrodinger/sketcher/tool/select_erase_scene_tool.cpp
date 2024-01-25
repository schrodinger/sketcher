#include "schrodinger/sketcher/tool/select_erase_scene_tool.h"

#include <QGraphicsItem>
#include <QList>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MolOps.h>

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/sgroup_item.h"
#include "schrodinger/sketcher/molviewer/non_molecular_item.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/subset.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * Create a QRectF with corners of the two given point.  Note that the
 * QRectF(QPointF, QPointF) constructor expects that the two arguments are the
 * upper-left and lower-right corners, respectively.  `rect_for_points` avoids
 * that requirement: any two corners in any order will result in a QRect with a
 * positive width and height.
 */
QRectF rect_for_points(const QPointF& a, const QPointF& b)
{
    qreal x = qMin(a.x(), b.x());
    qreal y = qMin(a.y(), b.y());
    qreal width = qAbs(a.x() - b.x());
    qreal height = qAbs(a.y() - b.y());
    return QRectF(x, y, width, height);
}

} // unnamed namespace

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
void SelectSceneTool<T>::onDragStart(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onDragStart(event);
    m_select_item.setVisible(true);
}

template <typename T>
void SelectSceneTool<T>::onDragMove(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onDragMove(event);
    auto items = m_scene->getCollidingItemsUsingBondMidpoints(&m_select_item);
    m_predictive_highlighting_item.highlightItems(items);
}

template <typename T>
void SelectSceneTool<T>::onDragRelease(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onDragRelease(event);
    m_select_item.setVisible(false);
    m_predictive_highlighting_item.clearHighlightingPath();
    auto items = m_scene->getCollidingItemsUsingBondMidpoints(&m_select_item);
    onSelectionMade(items, event);
}

template <typename T>
void SelectSceneTool<T>::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onMouseClick(event);
    QGraphicsItem* item = m_scene->getTopInteractiveItemAt(
        event->scenePos(), InteractiveItemFlag::ALL);
    if (item != nullptr) {
        onSelectionMade({item}, event);
    } else {
        onSelectionMade({}, event);
    }
}

template <typename T> void
SelectSceneTool<T>::onMouseDoubleClick(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onMouseDoubleClick(event);
    QGraphicsItem* item = m_scene->getTopInteractiveItemAt(
        event->scenePos(), InteractiveItemFlag::MOLECULAR);
    if (item == nullptr) {
        // double click wasn't on the molecule, so we do nothing
        return;
    } else if (auto* s_group_item = qgraphicsitem_cast<SGroupItem*>(item)) {
        // TODO: double-click on an s-group should open up the dialog
        return;
    }

    // double click was on the molecule, so select all connected atoms and bonds
    const RDKit::Atom* atom = nullptr;
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        atom = atom_item->getAtom();
    } else if (auto* bond_item = qgraphicsitem_cast<BondItem*>(item)) {
        atom = bond_item->getBond()->getBeginAtom();
    }
    auto [atoms_to_select, bonds_to_select] =
        get_connected_atoms_and_bonds(atom);
    auto select_mode = getSelectMode(event);
    if (select_mode == SelectMode::TOGGLE) {
        // Toggle behaves strangely because the first click gets processed
        // separately and selects the atom or bond, so we just treat a
        // Ctrl-click as a regular click here
        select_mode = SelectMode::SELECT_ONLY;
    }
    m_mol_model->select(atoms_to_select, bonds_to_select, {}, {}, select_mode);
}

template <typename T>
std::vector<QGraphicsItem*> SelectSceneTool<T>::getGraphicsItems()
{
    auto items = SceneToolWithPredictiveHighlighting::getGraphicsItems();
    items.push_back(&m_select_item);
    return items;
}

template <typename T> SelectMode
SelectSceneTool<T>::getSelectMode(QGraphicsSceneMouseEvent* const event) const
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
std::tuple<std::unordered_set<const RDKit::Atom*>,
           std::unordered_set<const RDKit::Bond*>,
           std::unordered_set<const RDKit::SubstanceGroup*>,

           std::unordered_set<const NonMolecularObject*>>
SelectSceneTool<T>::getModelObjectsForGraphicsItems(
    const QList<QGraphicsItem*>& items) const
{
    std::unordered_set<const RDKit::Atom*> atoms;
    std::unordered_set<const RDKit::Bond*> bonds;
    std::unordered_set<const RDKit::SubstanceGroup*> sgroups;
    std::unordered_set<const NonMolecularObject*> non_molecular_objects;
    for (auto cur_item : items) {
        if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(cur_item)) {
            atoms.insert(atom_item->getAtom());
        } else if (auto* bond_item = qgraphicsitem_cast<BondItem*>(cur_item)) {
            bonds.insert(bond_item->getBond());
        } else if (auto* sgroup_item =
                       qgraphicsitem_cast<SGroupItem*>(cur_item)) {
            sgroups.insert(sgroup_item->getSubstanceGroup());
        } else if (auto* non_molecular_item =
                       qgraphicsitem_cast<NonMolecularItem*>(cur_item)) {
            non_molecular_objects.insert(
                non_molecular_item->getNonMolecularObject());
        }
    }
    return {atoms, bonds, sgroups, non_molecular_objects};
}

template <typename T>
void SelectSceneTool<T>::onSelectionMade(const QList<QGraphicsItem*>& items,
                                         QGraphicsSceneMouseEvent* const event)
{
    auto select_mode = getSelectMode(event);
    auto [atoms, bonds, sgroups, non_molecular_objects] =
        getModelObjectsForGraphicsItems(items);
    m_mol_model->select(atoms, bonds, sgroups, non_molecular_objects,
                        select_mode);
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

QPixmap LassoSelectSceneTool::getCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/select_lasso.svg");
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

template <typename T> QPixmap ShapeSelectSceneTool<T>::getCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/select_square.svg");
}

EraseSceneTool::EraseSceneTool(Scene* scene, MolModel* mol_model) :
    RectSelectSceneTool(scene, mol_model)
{
}

void EraseSceneTool::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
    QGraphicsItem* item = m_scene->getTopInteractiveItemAt(
        event->scenePos(), InteractiveItemFlag::ALL);
    // if clicking a triple or double bond, do not erase but decrease bond order
    std::unordered_map<RDKit::Bond::BondType, BondTool> bond_tool_map = {
        {RDKit::Bond::TRIPLE, BondTool::DOUBLE},
        {RDKit::Bond::DOUBLE, BondTool::SINGLE},
    };
    if (BondItem* bond = qgraphicsitem_cast<BondItem*>(item)) {
        auto bond_order = bond->getBond()->getBondType();
        if (bond_tool_map.count(bond_order)) {
            m_mol_model->mutateBonds({bond->getBond()},
                                     bond_tool_map[bond_order]);
            // return early to avoid erasing the bond
            return;
        }
    }
    RectSelectSceneTool::onMouseClick(event);
}

void EraseSceneTool::onSelectionMade(const QList<QGraphicsItem*>& items,
                                     QGraphicsSceneMouseEvent* const event)
{
    // immediately clear the predictive highlighting, since the highlighted
    // items won't exist after the removeAtomsAndBonds call
    m_predictive_highlighting_item.clearHighlightingPath();
    auto [atoms, bonds, sgroups, non_molecular_objects] =
        getModelObjectsForGraphicsItems(items);
    m_mol_model->remove(atoms, bonds, sgroups, non_molecular_objects);
}

void EraseSceneTool::onMouseDoubleClick(QGraphicsSceneMouseEvent* const event)
{
    // disable the select tool double-click behavior
}

QPixmap EraseSceneTool::getCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/mode_erase.svg",
                                /* recolor = */ false);
}

} // namespace sketcher
} // namespace schrodinger
