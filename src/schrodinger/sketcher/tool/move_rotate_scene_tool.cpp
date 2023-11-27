#include "schrodinger/sketcher/tool/move_rotate_scene_tool.h"

#include <algorithm>

#include <GraphMol/MolOps.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"

namespace schrodinger
{
namespace sketcher
{

MoveSelectionItem::MoveSelectionItem() : QGraphicsRectItem()
{
    setPen(SELECTION_OUTLINE_COLOR);
    setZValue(static_cast<qreal>(ZOrder::SELECTION_HIGHLIGHTING));
}

MergeHintItem::MergeHintItem(QGraphicsItem* parent) : QGraphicsItemGroup(parent)
{
    setZValue(static_cast<qreal>(ZOrder::HINT));
    m_circle_pen.setColor(STRUCTURE_HINT_COLOR);
    m_circle_pen.setWidthF(DRAG_MERGE_HINT_WIDTH);
}

void MergeHintItem::setCoordinates(std::vector<QPointF> centers)
{
    clear();
    for (auto cur_center : centers) {
        auto circle = new QGraphicsEllipseItem(
            cur_center.x() - DRAG_MERGE_HINT_RADIUS,
            cur_center.y() - DRAG_MERGE_HINT_RADIUS, 2 * DRAG_MERGE_HINT_RADIUS,
            2 * DRAG_MERGE_HINT_RADIUS);
        circle->setPen(m_circle_pen);
        addToGroup(circle);
    }
}

void MergeHintItem::clear()
{
    for (auto* child : childItems()) {
        removeFromGroup(child);
        delete child;
    }
}

MoveRotateSceneTool::MoveRotateSceneTool(Scene* scene, MolModel* mol_model) :
    AbstractSceneTool(scene, mol_model)
{
    updateRotationItem();
    updateMoveSelectionItem();
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

    QString description;
    if (m_rotation_item.isInsideHandle(event->scenePos())) {
        m_action = Action::ROTATE;
        description = "Rotate";
        setObjectsToMove(selected_atoms, selected_non_mol_objs);
    } else {
        m_action = Action::TRANSLATE;
        description = "Translate";
        if (m_move_selection_item.rect().contains(event->scenePos())) {
            setObjectsToMove(selected_atoms, selected_non_mol_objs);
        }
    }
    // begin an undo macro in case we have to merge atoms at the end of the drag
    m_mol_model->beginUndoMacro(description);
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
    updateMergeHintItem();
    updateMoveSelectionItem();
}

void MoveRotateSceneTool::onDragRelease(QGraphicsSceneMouseEvent* event)
{
    mergeOverlappingAtoms();
    m_mol_model->endUndoMacro();
    m_action = Action::NONE;
    setObjectsToMove({}, {});
    updateMergeHintItem();
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
    return {&m_rotation_item, &m_move_selection_item, &m_merge_hint_item};
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
    updateRotationItem();
    updateMoveSelectionItem();
}

QPixmap MoveRotateSceneTool::getCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/select_move_rotate.svg");
}

std::vector<std::pair<unsigned int, unsigned int>>
MoveRotateSceneTool::getOverlappingAtomIdxs()
{
    return get_overlapping_atom_idxs(m_mol_model->getMol(), m_atoms_to_move);
}

void MoveRotateSceneTool::updateMergeHintItem()
{
    auto overlapping_idxs = getOverlappingAtomIdxs();
    std::vector<QPointF> circle_centers;
    const auto& conf = m_mol_model->getMol()->getConformer();
    // for each atom pair, get the coordinates of the stationary atom
    std::transform(overlapping_idxs.begin(), overlapping_idxs.end(),
                   std::back_inserter(circle_centers),
                   [&conf](std::pair<unsigned int, unsigned int> idxs) {
                       auto [to_move_idx, stationary_idx] = idxs;
                       const auto& stationary_coords =
                           conf.getAtomPos(stationary_idx);
                       return to_scene_xy(stationary_coords);
                   });
    m_merge_hint_item.setCoordinates(circle_centers);
}

void MoveRotateSceneTool::mergeOverlappingAtoms()
{
    auto overlapping_idxs = getOverlappingAtomIdxs();
    auto* mol = m_mol_model->getMol();
    std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>>
        overlapping_atoms;
    // convert from atom indices to Atom*
    std::transform(overlapping_idxs.begin(), overlapping_idxs.end(),
                   std::back_inserter(overlapping_atoms),
                   [&mol](std::pair<unsigned int, unsigned int> idxs) {
                       return std::make_pair(mol->getAtomWithIdx(idxs.first),
                                             mol->getAtomWithIdx(idxs.second));
                   });
    m_mol_model->mergeAtoms(overlapping_atoms);
}

std::vector<std::pair<unsigned int, unsigned int>> get_overlapping_atom_idxs(
    const RDKit::ROMol* const mol,
    const std::unordered_set<const RDKit::Atom*> atoms_to_move)
{
    std::vector<std::pair<unsigned int, unsigned int>> overlaps;
    if (atoms_to_move.empty()) {
        // nothing is being moved, so there are no overlaps to report
        return overlaps;
    }
    // figure out the indices of the moving and stationary atoms
    std::unordered_set<unsigned int> to_move_indices;
    std::unordered_set<unsigned int> stationary_indices;
    std::transform(atoms_to_move.begin(), atoms_to_move.end(),
                   std::inserter(to_move_indices, to_move_indices.begin()),
                   [](auto atom) { return atom->getIdx(); });
    for (unsigned int i = 0; i < mol->getNumAtoms(); ++i) {
        if (!to_move_indices.count(i)) {
            stationary_indices.insert(i);
        }
    }

    // figure out which atoms are actually overlapping
    std::vector<int> frag_map;
    const auto& conf = mol->getConformer();
    RDKit::MolOps::getMolFrags(*mol, frag_map);
    for (auto to_move_idx : to_move_indices) {
        auto to_move_frag_num = frag_map[to_move_idx];
        const auto& to_move_coords = conf.getAtomPos(to_move_idx);
        // we're looking for the closest overlapping atom
        int idx_to_merge = -1;
        double best_dist_sq = MAX_DIST_SQ_FOR_DRAG_MERGE;
        for (auto stationary_idx : stationary_indices) {
            if (frag_map[stationary_idx] == to_move_frag_num) {
                // these two atoms are from the same molecule, so we never
                // consider them overlapping
                continue;
            }
            const auto& stationary_coords = conf.getAtomPos(stationary_idx);
            auto dist_sq = (stationary_coords - to_move_coords).lengthSq();
            if (dist_sq <= best_dist_sq) {
                idx_to_merge = stationary_idx;
                best_dist_sq = dist_sq;
            }
        }
        if (idx_to_merge >= 0) {
            // there was an overlapping atom, so record the pair
            overlaps.emplace_back(to_move_idx, idx_to_merge);
        }
    }
    return overlaps;
}

} // namespace sketcher
} // namespace schrodinger