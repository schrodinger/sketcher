#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

#include <algorithm>

#include <GraphMol/MolOps.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{

AngleTextItem::AngleTextItem() : QGraphicsSimpleTextItem()
{
    setZValue(static_cast<qreal>(ZOrder::ROTATION_HANDLE));
    setBrush(ROTATION_ITEM_TEXT_COLOR);
    auto angle_font = font();
    angle_font.setPointSize(ROTATION_ITEM_FONT_SIZE);
    setFont(angle_font);
}

void AngleTextItem::centerOn(const QPointF& point)
{
    setPos(point - boundingRect().center());
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

StandardSceneToolBase::StandardSceneToolBase(Scene* scene,
                                             MolModel* mol_model) :
    AbstractSceneTool(scene, mol_model)
{
}

void StandardSceneToolBase::updateColorsAfterBackgroundColorChange(
    bool is_dark_mode)
{
    auto color = is_dark_mode ? PREDICTIVE_HIGHLIGHTING_COLOR_DARK_BG
                              : PREDICTIVE_HIGHLIGHTING_COLOR;
    m_predictive_highlighting_item.setPen(color);
    m_predictive_highlighting_item.setBrush(color);
    m_angle_text_item.setPen(is_dark_mode ? ANNOTATION_COLOR_DARK_BG
                                          : ANNOTATION_COLOR);
}

void StandardSceneToolBase::onRightButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    if (m_allow_context_menu) {
        emit contextMenuRequested(event);
    }
    AbstractSceneTool::onRightButtonClick(event);
}

void StandardSceneToolBase::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    if (m_mouse_pressed) {
        m_predictive_highlighting_item.clearHighlightingPath();
    } else if (m_highlight_types != InteractiveItemFlag::NONE) {
        QPointF scene_pos = event->scenePos();
        QGraphicsItem* item =
            m_scene->getTopInteractiveItemAt(scene_pos, m_highlight_types);
        m_predictive_highlighting_item.highlightItem(item);
    }
    AbstractSceneTool::onMouseMove(event);
}

void StandardSceneToolBase::onMiddleButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
    // if a selection is present, rotate only the selection
    setCurrentSelectionAsObjectsToMove();
    emit newCursorHintRequested(m_rotate_cursor_hint);
    emit atomDragStarted();
    m_mol_model->beginUndoMacro("Rotate");
    AbstractSceneTool::onMiddleButtonDragStart(event);
}

void StandardSceneToolBase::onMiddleButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    auto center_of_rotation = findPivotPointForRotation();
    rotate(event, to_scene_xy(center_of_rotation), m_atoms_to_move,
           m_non_mol_objs_to_move);
    AbstractSceneTool::onMiddleButtonDragMove(event);
}

void StandardSceneToolBase::onMiddleButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    finishDrag();
    m_angle_text_item.setText("");
    AbstractSceneTool::onMiddleButtonDragRelease(event);
}

void StandardSceneToolBase::onRightButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
    // if the translate drag starts inside a selection, move only the selection
    if (m_scene->getSelectionRect().contains(m_mouse_press_scene_pos)) {
        setCurrentSelectionAsObjectsToMove();
    }
    emit newCursorHintRequested(m_translate_cursor_hint);
    emit atomDragStarted();
    m_mol_model->beginUndoMacro("Translate");
    AbstractSceneTool::onRightButtonDragStart(event);
}

void StandardSceneToolBase::onRightButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    translate(event, m_atoms_to_move, m_non_mol_objs_to_move);
    updateMergeHintItem();
    AbstractSceneTool::onRightButtonDragMove(event);
}

void StandardSceneToolBase::onRightButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    finishDrag();
    updateMergeHintItem();
    AbstractSceneTool::onRightButtonDragRelease(event);
}

void StandardSceneToolBase::finishDrag()
{
    emit newCursorHintRequested(getDefaultCursorPixmap());
    mergeOverlappingAtoms();
    m_mol_model->endUndoMacro();
    setObjectsToMove({}, {});
}

void StandardSceneToolBase::setObjectsToMove(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const NonMolecularObject*>& non_mol_objects)
{
    m_atoms_to_move = atoms;
    m_non_mol_objs_to_move = non_mol_objects;
}

void StandardSceneToolBase::setCurrentSelectionAsObjectsToMove()
{
    auto selected_atoms = m_mol_model->getSelectedAtoms();
    auto selected_non_mol_objs = m_mol_model->getSelectedNonMolecularObjects();
    if (!selected_atoms.empty() || !selected_non_mol_objs.empty()) {
        setObjectsToMove(selected_atoms, selected_non_mol_objs);
    }
}

void StandardSceneToolBase::rotate(
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
    updateGraphicsItemsDuringRotation(event);
}

void StandardSceneToolBase::translate(
    QGraphicsSceneMouseEvent* const event,
    const std::unordered_set<const RDKit::Atom*>& atoms_to_move,
    const std::unordered_set<const NonMolecularObject*>& non_mol_objs_to_move)
{
    // if no objects are specified, move the viewport instead
    if (atoms_to_move.empty() && non_mol_objs_to_move.empty()) {
        // we can't use scenePos() here because the scene gets moved by the
        // viewport translation, so passing screen coordiates instead
        m_scene->viewportTranslationRequested(event->lastScreenPos(),
                                              event->screenPos());
    } else {
        auto distance = event->scenePos() - event->lastScenePos();
        m_mol_model->translateByVector(to_mol_xy(distance), atoms_to_move,
                                       non_mol_objs_to_move);
        updateGraphicsItemsDuringTranslation(distance);
    }
}

void StandardSceneToolBase::updateGraphicsItemsDuringRotation(
    QGraphicsSceneMouseEvent* const event)
{
    auto center = to_scene_xy(findPivotPointForRotation());
    auto angle_to_display =
        QLineF(center, event->scenePos())
            .angleTo(QLineF(center, m_drag_start_scene_pos));
    m_angle_text_item.setText(QString::number(angle_to_display, 'f', 0));
    m_angle_text_item.centerOn(center);
}

void StandardSceneToolBase::updateGraphicsItemsDuringTranslation(
    const QPointF& distance)
{
}

std::vector<QGraphicsItem*> StandardSceneToolBase::getGraphicsItems()
{
    return {&m_predictive_highlighting_item, &m_merge_hint_item,
            &m_angle_text_item};
}

RDGeom::Point3D StandardSceneToolBase::findPivotPointForRotation()
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

std::vector<std::pair<unsigned int, unsigned int>>
StandardSceneToolBase::getOverlappingAtomIdxs()
{
    return get_overlapping_atom_idxs(m_mol_model->getMol(), m_atoms_to_move);
}

void StandardSceneToolBase::updateMergeHintItem()
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

void StandardSceneToolBase::mergeOverlappingAtoms()
{
    auto overlapping_idxs = getOverlappingAtomIdxs();
    bool have_atoms_to_merge = !overlapping_idxs.empty();
    if (have_atoms_to_merge) {
        auto* mol = m_mol_model->getMol();
        std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>>
            overlapping_atoms;
        // convert from atom indices to Atom*
        std::transform(overlapping_idxs.begin(), overlapping_idxs.end(),
                       std::back_inserter(overlapping_atoms),
                       [&mol](std::pair<unsigned int, unsigned int> idxs) {
                           return std::make_pair(
                               mol->getAtomWithIdx(idxs.first),
                               mol->getAtomWithIdx(idxs.second));
                       });
        m_mol_model->mergeAtoms(overlapping_atoms);
    }
    emit atomDragFinished(have_atoms_to_merge);
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