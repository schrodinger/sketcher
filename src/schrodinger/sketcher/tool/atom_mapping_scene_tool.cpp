#include "schrodinger/sketcher/tool/atom_mapping_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/rdkit/molops.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{

ArrowHeadItem::ArrowHeadItem()
{
    setZValue(static_cast<qreal>(ZOrder::HINT));
    auto pen = QPen(HINT_COLOR);
    pen.setStyle(Qt::SolidLine);
    setPen(pen);
    setBrush(QBrush(HINT_COLOR));
    setVisible(false);
}

ArrowLineItem::ArrowLineItem()
{
    setZValue(static_cast<qreal>(ZOrder::HINT));
    auto pen = QPen(HINT_COLOR);
    pen.setStyle(Qt::DotLine);
    setPen(pen);
    setVisible(false);
}

AtomMappingSceneTool::AtomMappingSceneTool(const MappingAction& action,
                                           Scene* scene, MolModel* mol_model) :
    SceneToolWithPredictiveHighlighting(scene, mol_model),
    m_action(action)
{
    m_highlight_types = InteractiveItemFlag::ATOM_NOT_R_NOT_AP;
}

void AtomMappingSceneTool::updateArrowItems(const QLineF& line)
{
    m_arrow_line_item.setLine(line);

    // arrow head
    auto perpendicular_line = line.normalVector().unitVector();
    auto perpendicular = (perpendicular_line.p2() - perpendicular_line.p1()) *
                         MAPPING_ARROW_HEAD_HALF_WIDTH;
    auto parallel_line = line.unitVector();
    auto parallel =
        (parallel_line.p2() - parallel_line.p1()) * MAPPING_ARROW_HEAD_LENGTH;
    QPolygonF arrow_head;
    auto end = line.p2();
    arrow_head << end << end - parallel + perpendicular
               << end - parallel - perpendicular << end;
    m_arrow_head_item.setPolygon(arrow_head);
}

void AtomMappingSceneTool::onMousePress(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onMousePress(event);
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(scene_pos, m_highlight_types);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        m_pressed_atom_item = atom_item;
    }
}

void AtomMappingSceneTool::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onMouseClick(event);
    if (m_action != MappingAction::REMOVE) {
        return;
    }

    if (m_pressed_atom_item == nullptr) {
        return;
    }
    int mapping_n = m_pressed_atom_item->getAtom()->getAtomMapNum();
    if (mapping_n < 1) {
        return;
    }
    std::unordered_set<const RDKit::Atom*> atoms_to_remove;
    auto mol = m_mol_model->getMol();
    auto all_atoms = mol->atoms();
    copy_if(all_atoms.begin(), all_atoms.end(),
            std::inserter(atoms_to_remove, atoms_to_remove.begin()),
            [mapping_n](const RDKit::Atom* atom) {
                return atom->getAtomMapNum() == mapping_n;
            });
    m_mol_model->setAtomMapping(atoms_to_remove, 0);
}

void AtomMappingSceneTool::onDragMove(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onDragMove(event);
    if (m_action != MappingAction::ADD) {
        return;
    }
    if (m_pressed_atom_item == nullptr) {
        return;
    }

    QPointF scene_pos = event->scenePos();
    QGraphicsItem* item =
        m_scene->getTopInteractiveItemAt(scene_pos, m_highlight_types);
    auto atom_item = qgraphicsitem_cast<AtomItem*>(item);

    // TODO:  in reactions avoid mapping within products or within reactants

    if (atom_item != m_pressed_atom_item) {
        m_release_atom_item = atom_item;
    }

    // snap the line to hovered atom if valid
    if (atom_item != nullptr) {
        scene_pos = atom_item->scenePos();
    }

    // shorten line to avoid overlapping with atoms
    QLineF line(m_pressed_atom_item->scenePos(), scene_pos);

    auto trim_line_to_atom = [](QLineF& line, AtomItem* atom_item) {
        auto subrects = atom_item->getSubrects();
        if (subrects.empty()) {
            // if no subrects (implicit Carbon), use a small one to avoid being
            // too close to the bond lines). trim_line_to_rect will add a 4
            // pixel margin around this
            auto SMALL = 0.1;
            subrects.push_back(QRectF(-SMALL, -SMALL, SMALL * 2, SMALL * 2));
        }
        for (const QRectF& subrect : subrects) {

            trim_line_to_rect(line, subrect.translated(atom_item->scenePos()));
        }
    };

    trim_line_to_atom(line, m_pressed_atom_item);
    if (atom_item != nullptr) {
        trim_line_to_atom(line, atom_item);
    }
    updateArrowItems(line);
    m_arrow_head_item.setVisible(true);
    m_arrow_line_item.setVisible(true);
}

void AtomMappingSceneTool::onDragRelease(QGraphicsSceneMouseEvent* const event)
{
    SceneToolWithPredictiveHighlighting::onDragRelease(event);
    if (m_pressed_atom_item && m_release_atom_item) {
        // if the clicked atom is mapped, use that mapping for both atoms
        int mapping = m_pressed_atom_item->getAtom()->getAtomMapNum();
        if (mapping < 1) {
            mapping = findLowestAvailableMappingNumber();
        }
        m_mol_model->setAtomMapping(
            {m_pressed_atom_item->getAtom(), m_release_atom_item->getAtom()},
            mapping);
    }
    m_pressed_atom_item = nullptr;
    m_release_atom_item = nullptr;
    m_arrow_head_item.setVisible(false);
    m_arrow_line_item.setVisible(false);
}

std::vector<QGraphicsItem*> AtomMappingSceneTool::getGraphicsItems()
{
    auto items = SceneToolWithPredictiveHighlighting::getGraphicsItems();
    items.push_back(&m_arrow_line_item);
    items.push_back(&m_arrow_head_item);

    return items;
}

int AtomMappingSceneTool::findLowestAvailableMappingNumber() const
{
    auto mol = m_mol_model->getMol();
    auto all_atoms = mol->atoms();
    std::set<int> used_mapping_numbers;
    for (auto atom : all_atoms) {
        used_mapping_numbers.insert(atom->getAtomMapNum());
    }
    int lowest_available_mapping_number = 1;

    while (std::find(used_mapping_numbers.begin(), used_mapping_numbers.end(),
                     lowest_available_mapping_number) !=
           used_mapping_numbers.end()) {
        lowest_available_mapping_number++;
    }
    return lowest_available_mapping_number;
}

} // namespace sketcher
} // namespace schrodinger
