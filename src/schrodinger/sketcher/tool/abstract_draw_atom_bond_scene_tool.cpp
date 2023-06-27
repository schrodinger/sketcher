#include "schrodinger/sketcher/tool/abstract_draw_atom_bond_scene_tool.h"

#include <cmath>

#include <QtMath>
#include <QPen>

#include <Geometry/point.h>
#include <GraphMol/ROMol.h>

#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"

namespace schrodinger
{
namespace sketcher
{

HintBondItem::HintBondItem(QGraphicsItem* parent) : QGraphicsLineItem(parent)
{
    setZValue(static_cast<qreal>(ZOrder::HINT));
    setPen(QPen(HINT_COLOR));
    setVisible(false);
}

AbstractDrawSceneTool::AbstractDrawSceneTool(Scene* scene,
                                             MolModel* mol_model) :
    SceneToolWithPredictiveHighlighting(scene, mol_model)
{
    m_highlight_attachment_points = false;
}

std::vector<QGraphicsItem*> AbstractDrawSceneTool::getGraphicsItems()
{
    auto items = SceneToolWithPredictiveHighlighting::getGraphicsItems();
    items.push_back(&m_hint_bond_item);
    return items;
}

void AbstractDrawSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag logic is handled in onDragMove
        return;
    }
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(
        scene_pos, /*include_attachment_points = */ false);
    bool drew_hint = false;
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        const RDKit::Atom* atom = atom_item->getAtom();
        if (shouldDrawBondForClickOnAtom(atom)) {
            auto [mol_pos, scene_pos, end_atom] = getDefaultBondPosition(atom);
            updateHintBondPath(atom_item->pos(), scene_pos);
            drew_hint = true;
        }
    }
    setHintBondVisible(drew_hint);
}

void AbstractDrawSceneTool::onMouseRelease(
    QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onMouseRelease(event);
    setHintBondVisible(false);
}

void AbstractDrawSceneTool::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onMouseClick(event);
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(
        scene_pos, /*include_attachment_points = */ false);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        const RDKit::Atom* atom = atom_item->getAtom();
        onAtomClicked(atom);
    } else if (auto* bond_item = qgraphicsitem_cast<BondItem*>(item)) {
        const RDKit::Bond* bond = bond_item->getBond();
        onBondClicked(bond);
    } else {
        RDGeom::Point3D mol_pos = to_mol_xy(scene_pos);
        onEmptySpaceClicked(mol_pos);
    }
}

void AbstractDrawSceneTool::onDragStart(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onDragStart(event);
    auto [should_drag, start_pos, start_atom] = getDragStartInfo();
    setHintBondVisible(should_drag);
}

void AbstractDrawSceneTool::onDragMove(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onDragMove(event);
    auto [should_drag, start_pos, start_atom] = getDragStartInfo();
    if (should_drag) {
        auto [end_pos, end_atom] =
            getBondEndInMousedDirection(start_pos, event->scenePos());
        updateHintBondPath(start_pos, end_pos);
    }
}

void AbstractDrawSceneTool::onDragRelease(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onDragRelease(event);
    auto [should_drag, start_pos, start_atom] = getDragStartInfo();
    if (should_drag) {
        auto [end_pos, end_atom] =
            getBondEndInMousedDirection(start_pos, event->scenePos());
        if (start_atom && end_atom) {
            addBond(start_atom, end_atom);
        } else if (start_atom) {
            addAtom(to_mol_xy(end_pos), start_atom);
        } else if (end_atom) {
            addAtom(to_mol_xy(start_pos), end_atom);
        } else {
            addTwoBoundAtoms(to_mol_xy(start_pos), to_mol_xy(end_pos));
        }
    }
}

void AbstractDrawSceneTool::setHintBondVisible(bool visible)
{
    m_hint_bond_item.setVisible(visible);
}

void AbstractDrawSceneTool::updateHintBondPath(const QPointF& start,
                                               const QPointF& end)
{
    QLineF line = QLineF(start, end);
    m_hint_bond_item.setLine(line);
}

void AbstractDrawSceneTool::cycleBond(const RDKit::Bond* const bond)
{
    auto bond_type = bond->getBondType();
    if (bond_type == RDKit::Bond::BondType::SINGLE) {
        m_mol_model->mutateBond(bond, RDKit::Bond::BondType::DOUBLE);
    } else if (bond_type == RDKit::Bond::BondType::DOUBLE) {
        m_mol_model->mutateBond(bond, RDKit::Bond::BondType::TRIPLE);
    } else {
        m_mol_model->mutateBond(bond, RDKit::Bond::BondType::SINGLE);
    }
}

std::pair<RDGeom::Point3D, QPointF>
AbstractDrawSceneTool::getInitialDefaultBondPosition(
    const RDKit::Atom* const atom, const float bond_length) const
{
    auto positions = get_relative_positions_of_atom_neighbors(atom);
    RDGeom::Point3D best_offset = best_placing_around_origin(positions);
    best_offset *= bond_length;
    auto& conf = m_mol_model->getMol()->getConformer();
    RDGeom::Point3D atom_pos = conf.getAtomPos(atom->getIdx());
    RDGeom::Point3D bond_end = atom_pos + best_offset;
    QPointF scene_bond_end = to_scene_xy(bond_end);
    return {bond_end, scene_bond_end};
}

std::tuple<RDGeom::Point3D, QPointF, const RDKit::Atom*>
AbstractDrawSceneTool::getDefaultBondPosition(
    const RDKit::Atom* const atom) const
{
    auto [bond_end, scene_bond_end] = getInitialDefaultBondPosition(atom);
    AbstractGraphicsItem* item = m_scene->getTopInteractiveItemAt(
        scene_bond_end, /*include_attachment_points = */ false);
    const RDKit::Atom* atom_at_bond_end = nullptr;
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        atom_at_bond_end = atom_item->getAtom();
        scene_bond_end = atom_item->pos();
    }
    return {bond_end, scene_bond_end, atom_at_bond_end};
}

std::tuple<bool, QPointF, const RDKit::Atom*>
AbstractDrawSceneTool::getDragStartInfo() const
{
    auto* item = m_scene->getTopInteractiveItemAt(
        m_mouse_press_scene_pos, /*include_attachment_points = */ false);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        return {true, atom_item->pos(), atom_item->getAtom()};
    } else if (qgraphicsitem_cast<BondItem*>(item)) {
        return {false, QPointF(), nullptr};
    } else {
        return {true, m_mouse_press_scene_pos, nullptr};
    }
}

std::pair<QPointF, const RDKit::Atom*>
AbstractDrawSceneTool::getBondEndInMousedDirection(
    const QPointF& start, const QPointF& mouse_pos) const
{
    auto* item = m_scene->getTopInteractiveItemAt(
        mouse_pos, /*include_attachment_points = */ false);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        return {atom_item->pos(), atom_item->getAtom()};
    }
    QPointF bond_offset =
        getDefaultBondOffsetInMousedDirection(start, mouse_pos);
    QPointF bond_end = start + bond_offset;
    item = m_scene->getTopInteractiveItemAt(
        bond_end,
        /*include_attachment_points = */ false);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        return {atom_item->pos(), atom_item->getAtom()};
    }
    return {bond_end, nullptr};
}

QPointF AbstractDrawSceneTool::getDefaultBondOffsetInMousedDirection(
    const QPointF& start, const QPointF& mouse_pos) const
{
    qreal angle = QLineF(start, mouse_pos).angle();
    angle = qDegreesToRadians(angle);
    int rounded = std::round(angle * 6.0 / M_PI);
    qreal new_angle = rounded / 6.0 * M_PI;
    return QPointF(VIEW_SCALE * qCos(new_angle), -VIEW_SCALE * qSin(new_angle));
}

void AbstractDrawSceneTool::onAtomClicked(const RDKit::Atom* const atom)
{
    if (shouldDrawBondForClickOnAtom(atom)) {
        auto [mol_pos, scene_pos, end_atom] = getDefaultBondPosition(atom);
        if (end_atom) {
            // there's already an atom where the default bond position ends
            addBond(atom, end_atom);
        } else {
            addAtom(mol_pos, atom);
        }
    } else {
        mutateAtom(atom);
    }
}

void AbstractDrawSceneTool::mutateAtom(const RDKit::Atom* const atom)
{
    // This method is intentionally left blank
}

} // namespace sketcher
} // namespace schrodinger
