#include "schrodinger/sketcher/tool/abstract_draw_atom_bond_scene_tool.h"

#include <cmath>

#include <QtMath>
#include <QPen>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>

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
    setPen(QPen(STRUCTURE_HINT_COLOR));
    setVisible(false);
}

AbstractDrawSceneTool::AbstractDrawSceneTool(Scene* scene,
                                             MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model)
{
    m_highlight_types = InteractiveItemFlag::MOLECULAR_NOT_AP;
}

std::vector<QGraphicsItem*> AbstractDrawSceneTool::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
    items.push_back(&m_hint_bond_item);
    return items;
}

void AbstractDrawSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag logic is handled in onDragMove
        return;
    }
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(
        scene_pos, InteractiveItemFlag::MOLECULAR_NOT_AP);
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

void AbstractDrawSceneTool::onLeftButtonRelease(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonRelease(event);
    setHintBondVisible(false);
}

void AbstractDrawSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonClick(event);
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(
        scene_pos, InteractiveItemFlag::MOLECULAR_NOT_AP);
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

void AbstractDrawSceneTool::onLeftButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragStart(event);
    auto [should_drag, start_pos, start_atom] = getDragStartInfo();
    setHintBondVisible(should_drag);
}

void AbstractDrawSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragMove(event);
    auto [should_drag, start_pos, start_atom] = getDragStartInfo();
    if (should_drag) {
        auto [end_pos, end_atom] = getBondEndInMousedDirection(
            start_pos, start_atom, event->scenePos());
        updateHintBondPath(start_pos, end_pos);
    }
}

void AbstractDrawSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragRelease(event);
    auto [should_drag, start_pos, start_atom] = getDragStartInfo();
    if (should_drag) {
        auto [end_pos, end_atom] = getBondEndInMousedDirection(
            start_pos, start_atom, event->scenePos());
        if (start_atom && end_atom) {
            auto* mol = m_mol_model->getMol();
            auto* existing_bond = mol->getBondBetweenAtoms(start_atom->getIdx(),
                                                           end_atom->getIdx());
            if (existing_bond == nullptr) {
                addBond(start_atom, end_atom);
            } else {
                // The user dragged along an existing bond, so we respond as if
                // that bond was clicked
                onBondClicked(existing_bond);
            }
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
        m_mol_model->mutateBonds({bond}, BondTool::DOUBLE);
    } else if (bond_type == RDKit::Bond::BondType::DOUBLE) {
        m_mol_model->mutateBonds({bond}, BondTool::TRIPLE);
    } else {
        m_mol_model->mutateBonds({bond}, BondTool::SINGLE);
    }
}

std::pair<RDGeom::Point3D, QPointF>
AbstractDrawSceneTool::getInitialDefaultBondPosition(
    const RDKit::Atom* const atom) const
{
    auto positions = get_relative_positions_of_atom_neighbors(atom);
    RDGeom::Point3D best_offset = best_placing_around_origin(positions);
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
        scene_bond_end, InteractiveItemFlag::MOLECULAR_NOT_AP);
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
        m_mouse_press_scene_pos, InteractiveItemFlag::MOLECULAR_NOT_AP);
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
    const QPointF& start_pos, const RDKit::Atom* const start_atom,
    const QPointF& mouse_pos) const
{
    auto* item = m_scene->getTopInteractiveItemAt(
        mouse_pos, InteractiveItemFlag::ATOM_NOT_AP);
    bool force_default_direction = false;
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        auto* atom = atom_item->getAtom();
        if (start_atom != nullptr && atom == start_atom) {
            // The user is mousing over the start atom, so we'll use the default
            // bond direction below
            force_default_direction = true;
        } else {
            // The user is mousing over an atom, but not the start atom
            return {atom_item->pos(), atom};
        }
    }

    QPointF bond_end;
    if (force_default_direction) {
        // The user is mousing over the start atom
        bond_end = getInitialDefaultBondPosition(start_atom).second;
    } else {
        QPointF bond_offset =
            getDefaultBondOffsetInMousedDirection(start_pos, mouse_pos);
        bond_end = start_pos + bond_offset;
    }
    item = m_scene->getTopInteractiveItemAt(bond_end,
                                            InteractiveItemFlag::ATOM_NOT_AP);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        // there's already an atom where we're planning on ending the bond
        return {atom_item->pos(), atom_item->getAtom()};
    }
    return {bond_end, nullptr};
}

QPointF AbstractDrawSceneTool::getDefaultBondOffsetInMousedDirection(
    const QPointF& start, const QPointF& mouse_pos) const
{
    qreal new_angle = get_rounded_angle_radians(start, mouse_pos);
    return to_scene_xy(RDGeom::Point3D(BOND_LENGTH * qCos(new_angle),
                                       BOND_LENGTH * qSin(new_angle), 0));
}

void AbstractDrawSceneTool::onAtomClicked(const RDKit::Atom* const atom)
{
    if (shouldDrawBondForClickOnAtom(atom)) {
        auto [mol_pos, scene_pos, end_atom] = getDefaultBondPosition(atom);
        if (end_atom) {
            // there's already an atom where the default bond position ends
            auto* mol = m_mol_model->getMol();
            auto* existing_bond =
                mol->getBondBetweenAtoms(atom->getIdx(), end_atom->getIdx());
            // If there's already an existing bond between these two atoms, then
            // it means that the start atom already has so many bonds that we
            // can't find enough free space to draw a new one.  If that's
            // happened, just ignore the click.  This normally doesn't happen
            // until the 25th bond, so there's no point in trying to find
            // reasonable chemistry.
            if (existing_bond == nullptr) {
                addBond(atom, end_atom);
            }
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
