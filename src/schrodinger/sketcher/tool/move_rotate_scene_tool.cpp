#include "schrodinger/sketcher/tool/move_rotate_scene_tool.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{

MoveSelectionItem::MoveSelectionItem() : QGraphicsRectItem()
{
    setPen(SELECTION_OUTLINE_COLOR);
    setZValue(static_cast<qreal>(ZOrder::SELECTION_HIGHLIGHTING));
}

MoveRotateSceneTool::MoveRotateSceneTool(Scene* scene, MolModel* mol_model) :
    AbstractSceneTool(scene, mol_model)
{
}

// TODO: remove this, for testing only until context menus are implemented
void MoveRotateSceneTool::onMouseClick(QGraphicsSceneMouseEvent* event)
{
    auto mol = m_mol_model->getMol();
    auto selected_atoms = m_mol_model->getSelectedAtoms();
    if (selected_atoms.empty()) {
        return;
    }
    // Find all bonds with one selected atom and one unselected atom
    std::unordered_set<const RDKit::Bond*> crossing_bonds;
    for (auto bond : mol->bonds()) {
        if (selected_atoms.count(bond->getBeginAtom()) !=
            selected_atoms.count(bond->getEndAtom())) {
            crossing_bonds.insert(bond);
        }
    }
    // If there is only one, use it as the point
    if (crossing_bonds.size() == 1) {
        auto bond = *crossing_bonds.begin();
        auto start_coord =
            mol->getConformer().getAtomPos(bond->getBeginAtom()->getIdx());
        auto end_coord =
            mol->getConformer().getAtomPos(bond->getEndAtom()->getIdx());
        m_mol_model->flipAroundSegment(start_coord, end_coord, selected_atoms);
    }
}

void MoveRotateSceneTool::onDragStart(QGraphicsSceneMouseEvent* event)
{
    auto selected_atoms = m_mol_model->getSelectedAtoms();
    auto selected_non_mol_objs = m_mol_model->getSelectedNonMolecularObjects();

    if (m_rotation_item.isInsideHandle(event->scenePos())) {
        m_action = Action::ROTATE;
        setObjectsToMove(selected_atoms, selected_non_mol_objs);
    } else {
        m_action = Action::TRANSLATE;
        if (m_move_selection_item.rect().contains(event->scenePos())) {
            setObjectsToMove(selected_atoms, selected_non_mol_objs);
        }
    }
}

void MoveRotateSceneTool::onDragMove(QGraphicsSceneMouseEvent* event)
{
    switch (m_action) {
        case Action::ROTATE:
            rotate(event, m_rotation_item.getPivotPoint(), m_atoms_to_move,
                   m_non_mol_objs_to_move);
            break;
        case Action::TRANSLATE:
            translate(event, m_atoms_to_move, m_non_mol_objs_to_move);
            break;
        default:
            break;
    }
    updateMoveSelectionItem();
}

void MoveRotateSceneTool::onDragRelease(QGraphicsSceneMouseEvent* event)
{
    m_action = Action::NONE;
    setObjectsToMove({}, {});
}

void MoveRotateSceneTool::setObjectsToMove(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const NonMolecularObject*>& non_mol_objects)
{
    m_atoms_to_move = atoms;
    m_non_mol_objs_to_move = non_mol_objects;
}

void MoveRotateSceneTool::rotate(
    QGraphicsSceneMouseEvent* const event, QPointF pivot_point,
    const std::unordered_set<const RDKit::Atom*>& atoms_to_move,
    const std::unordered_set<const NonMolecularObject*>& non_mol_objs_to_move)
{
    // calculate angle increment of the event,
    // considering the centroid of the molecule as the origin
    QLineF current_position_line(pivot_point, event->scenePos());
    QLineF last_position_line(pivot_point, event->lastScenePos());
    auto angle = last_position_line.angleTo(current_position_line);
    m_mol_model->rotateByAngle(angle, to_mol_xy(pivot_point), atoms_to_move,
                               non_mol_objs_to_move);
    rotateRotationItem(event);
}

void MoveRotateSceneTool::rotateRotationItem(QGraphicsSceneMouseEvent* event)
{
    auto center = m_rotation_item.getPivotPoint();
    auto prev_angle = QLineF(center, event->lastScenePos()).angle();
    auto angle = QLineF(center, event->scenePos()).angle();
    m_rotation_item.setArmAngle(m_rotation_item.getArmAngle() + angle -
                                prev_angle);
}

void MoveRotateSceneTool::translate(
    QGraphicsSceneMouseEvent* const event,
    const std::unordered_set<const RDKit::Atom*>& atoms_to_move,
    const std::unordered_set<const NonMolecularObject*>& non_mol_objs_to_move)
{
    // if no atoms are specified, move the viewport instead
    if (atoms_to_move.empty()) {
        // we can't use scenePos() here because the scene gets moved by the
        // viewport translation, so passing screen coordiates instead
        m_scene->viewportTranslationRequested(event->lastScreenPos(),
                                              event->screenPos());
    } else {
        auto distance = event->scenePos() - event->lastScenePos();
        m_mol_model->translateByVector(to_mol_xy(distance), atoms_to_move);
        m_rotation_item.setPivotPoint(m_rotation_item.getPivotPoint() +
                                      distance);
    }
}

void MoveRotateSceneTool::updateRotationItem()
{
    if (m_mol_model->isEmpty()) {
        return;
    }
    m_rotation_item.setArmAngle(0);
    auto pivot_point = findPivotPointForRotation();
    m_rotation_item.setPivotPoint(to_scene_xy(pivot_point));

    // Show the rotation handle if there is more than one selected atom or
    // non-molecular object, or if there is no selection and the molecule has
    // more than one atom or non-molecular object
    auto sel_atoms = m_mol_model->getSelectedAtoms();
    auto sel_nmo = m_mol_model->getSelectedNonMolecularObjects();
    int num_atoms = sel_atoms.empty() && sel_nmo.empty()
                        ? m_mol_model->getMol()->getNumAtoms() +
                              m_mol_model->getNonMolecularObjects().size()
                        : sel_atoms.size() + sel_nmo.size();
    m_rotation_item.setDrawRotationHandle(num_atoms > 1);
}

void MoveRotateSceneTool::updateMoveSelectionItem()
{
    m_move_selection_item.setRect(m_scene->getSelectionRect());
}

std::vector<QGraphicsItem*> MoveRotateSceneTool::getGraphicsItems()
{
    updateRotationItem();
    updateMoveSelectionItem();

    return {&m_rotation_item, &m_move_selection_item};
}

RDGeom::Point3D MoveRotateSceneTool::findPivotPointForRotation()
{
    // if there's no coordinates, return the origin
    if (m_mol_model->isEmpty()) {
        return RDGeom::Point3D(0, 0, 0);
    }

    auto mol = m_mol_model->getMol();
    auto selected_atoms = m_mol_model->getSelectedAtoms();
    auto selected_non_mol_objs = m_mol_model->getSelectedNonMolecularObjects();
    if (selected_atoms.empty() && selected_non_mol_objs.empty()) {
        return find_centroid(*mol, m_mol_model->getNonMolecularObjects());
    }

    // Find all bonds with one selected atom and one unselected atom
    std::unordered_set<const RDKit::Bond*> crossing_bonds;
    for (auto bond : mol->bonds()) {
        if (selected_atoms.count(bond->getBeginAtom()) !=
            selected_atoms.count(bond->getEndAtom())) {
            crossing_bonds.insert(bond);
        }
    }
    // If there is only one, use the selected atom as the pivot point
    if (crossing_bonds.size() == 1) {
        auto bond = *crossing_bonds.begin();
        auto begin_selected = selected_atoms.count(bond->getBeginAtom()) == 1;
        auto pivot_point_idx = begin_selected ? bond->getBeginAtom()->getIdx()
                                              : bond->getEndAtom()->getIdx();
        return mol->getConformer().getAtomPos(pivot_point_idx);
    }
    return find_centroid(*mol, selected_atoms, selected_non_mol_objs);
}

void MoveRotateSceneTool::onStructureUpdated()
{
    updateMoveSelectionItem();
}

} // namespace sketcher
} // namespace schrodinger