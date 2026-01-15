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
#include "schrodinger/sketcher/rdkit/monomer_connectors.h"
#include "schrodinger/sketcher/rdkit/subset.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

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
    StandardSceneToolBase(scene, mol_model)
{
    m_highlight_types = InteractiveItemFlag::ALL;
    m_select_item.setVisible(false);
}

template <typename T> void
SelectSceneTool<T>::updateColorsAfterBackgroundColorChange(bool is_dark_mode)
{
    m_select_item.setPen(is_dark_mode ? SELECT_TOOL_LINE_COLOR_DARK_BG
                                      : SELECT_TOOL_LINE_COLOR);
    StandardSceneToolBase::updateColorsAfterBackgroundColorChange(is_dark_mode);
}

template <typename T> void
SelectSceneTool<T>::onLeftButtonDragStart(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragStart(event);
    m_select_item.setVisible(true);
}

template <typename T> void
SelectSceneTool<T>::onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragMove(event);
    auto items = m_scene->getCollidingItemsUsingBondMidpoints(&m_select_item);
    m_predictive_highlighting_item.highlightItems(items);
}

template <typename T> void SelectSceneTool<T>::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragRelease(event);
    m_select_item.setVisible(false);
    m_predictive_highlighting_item.clearHighlightingPath();
    auto items = m_scene->getCollidingItemsUsingBondMidpoints(&m_select_item);
    onSelectionMade(items, event);
}

template <typename T> void
SelectSceneTool<T>::onLeftButtonClick(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonClick(event);
    QGraphicsItem* item = m_scene->getTopInteractiveItemAt(
        event->scenePos(), InteractiveItemFlag::ALL);
    if (item != nullptr) {
        onSelectionMade({item}, event);
    } else {
        onSelectionMade({}, event);
    }
}

template <typename T> void SelectSceneTool<T>::onLeftButtonDoubleClick(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDoubleClick(event);
    QGraphicsItem* item = m_scene->getTopInteractiveItemAt(
        event->scenePos(), InteractiveItemFlag::MOLECULAR_OR_MONOMERIC);
    if (item == nullptr) {
        // double click wasn't on the molecule, so we do nothing
        return;
    } else if (qgraphicsitem_cast<SGroupItem*>(item)) {
        // TODO: double-click on an s-group should open up the dialog
        return;
    }

    // double click was on the molecule, so select all connected atoms and bonds
    const RDKit::Atom* atom = nullptr;
    if (auto* atom_item = dynamic_cast<AbstractAtomOrMonomerItem*>(item)) {
        atom = atom_item->getAtom();
    } else if (auto* bond_item =
                   dynamic_cast<AbstractBondOrConnectorItem*>(item)) {
        atom = bond_item->getBond()->getBeginAtom();
    }
    auto [atoms_to_select, bonds_to_select] =
        get_connected_atoms_and_bonds(atom);
    std::unordered_set<const RDKit::Bond*> secondary_connections_to_select;
    // if any of the bonds have a secondary connection, select that as well
    std::copy_if(bonds_to_select.begin(), bonds_to_select.end(),
                 std::inserter(secondary_connections_to_select,
                               secondary_connections_to_select.end()),
                 contains_two_monomer_linkages);
    auto select_mode = getSelectMode(event);
    if (select_mode == SelectMode::TOGGLE) {
        // Toggle behaves strangely because the first click gets processed
        // separately and selects the atom or bond, so we just treat a
        // Ctrl-click as a regular click here
        select_mode = SelectMode::SELECT_ONLY;
    }
    m_mol_model->select(atoms_to_select, bonds_to_select,
                        secondary_connections_to_select, {}, {}, select_mode);
}

template <typename T> void
SelectSceneTool<T>::onRightButtonClick(QGraphicsSceneMouseEvent* const event)
{
    // when right-clicking on an atom or bond that's not part of the selection,
    // select it
    auto pos = event->scenePos();
    auto item = m_scene->getTopInteractiveItemAt(
        pos, InteractiveItemFlag::MOLECULAR_OR_MONOMERIC);
    if (item && !item->isSelected()) {
        auto [atoms, bonds, secondary_connections, sgroups,
              non_molecular_objects] =
            m_scene->getModelObjects(SceneSubset::HOVERED, &pos);
        m_mol_model->select(atoms, bonds, secondary_connections, sgroups,
                            non_molecular_objects, SelectMode::SELECT_ONLY);
    }
    StandardSceneToolBase::onRightButtonClick(event);
}

template <typename T>
std::vector<QGraphicsItem*> SelectSceneTool<T>::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
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
void SelectSceneTool<T>::onSelectionMade(const QList<QGraphicsItem*>& items,
                                         QGraphicsSceneMouseEvent* const event)
{
    auto select_mode = getSelectMode(event);
    auto [atoms, bonds, secondary_connections, sgroups, non_molecular_objects] =
        get_model_objects_for_graphics_items(items);
    m_mol_model->select(atoms, bonds, secondary_connections, sgroups,
                        non_molecular_objects, select_mode);
}

LassoSelectSceneTool::LassoSelectSceneTool(Scene* scene, MolModel* mol_model) :
    SelectSceneTool(scene, mol_model)
{
}

void LassoSelectSceneTool::onLeftButtonPress(QGraphicsSceneMouseEvent* event)
{
    SelectSceneTool<LassoSelectionItem>::onLeftButtonPress(event);
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

void LassoSelectSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* event)
{
    SelectSceneTool<LassoSelectionItem>::onLeftButtonDragRelease(event);
    m_select_item.clearPath();
}

QPixmap LassoSelectSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/select_lasso.svg");
}

template <typename T>
ShapeSelectSceneTool<T>::ShapeSelectSceneTool(Scene* scene,
                                              MolModel* mol_model) :
    SelectSceneTool<T>(scene, mol_model)
{
}

template <typename T> void
ShapeSelectSceneTool<T>::onLeftButtonDragMove(QGraphicsSceneMouseEvent* event)
{
    QRectF rect =
        QRectF(this->m_mouse_press_scene_pos, event->scenePos()).normalized();
    this->m_select_item.setRect(rect);
    SelectSceneTool<T>::onLeftButtonDragMove(event);
}

template <> QPixmap RectSelectSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/select_square.svg");
}

template <> QPixmap EllipseSelectSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/select_ellipse.svg");
}

EraseSceneTool::EraseSceneTool(Scene* scene, MolModel* mol_model) :
    RectSelectSceneTool(scene, mol_model)
{
}

void EraseSceneTool::onLeftButtonClick(QGraphicsSceneMouseEvent* const event)
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
    RectSelectSceneTool::onLeftButtonClick(event);
}

void EraseSceneTool::onSelectionMade(const QList<QGraphicsItem*>& items,
                                     QGraphicsSceneMouseEvent* const event)
{
    // immediately clear the predictive highlighting, since the highlighted
    // items won't exist after the removeAtomsAndBonds call
    m_predictive_highlighting_item.clearHighlightingPath();
    auto [atoms, bonds, secondary_connections, sgroups, non_molecular_objects] =
        get_model_objects_for_graphics_items(items);
    m_mol_model->remove(atoms, bonds, secondary_connections, sgroups,
                        non_molecular_objects);
}

void EraseSceneTool::onLeftButtonDoubleClick(
    QGraphicsSceneMouseEvent* const event)
{
    // disable the select tool double-click behavior
}

void EraseSceneTool::onRightButtonClick(QGraphicsSceneMouseEvent* const event)
{
    // disable the select tool right-click behavior
    StandardSceneToolBase::onRightButtonClick(event);
}

QPixmap EraseSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/mode_erase.svg",
                                /* recolor = */ false);
}

} // namespace sketcher
} // namespace schrodinger
