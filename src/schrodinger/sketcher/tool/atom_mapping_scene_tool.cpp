#include "schrodinger/sketcher/tool/atom_mapping_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"

namespace schrodinger
{
namespace sketcher
{

ArrowHeadItem::ArrowHeadItem()
{
    setZValue(static_cast<qreal>(ZOrder::HINT));
    auto pen = QPen(STRUCTURE_HINT_COLOR);
    pen.setStyle(Qt::SolidLine);
    setPen(pen);
    setBrush(QBrush(STRUCTURE_HINT_COLOR));
    setVisible(false);
}

ArrowLineItem::ArrowLineItem()
{
    setZValue(static_cast<qreal>(ZOrder::HINT));
    auto pen = QPen(STRUCTURE_HINT_COLOR);
    pen.setStyle(Qt::DotLine);
    setPen(pen);
    setVisible(false);
}

AtomMappingSceneTool::AtomMappingSceneTool(const MappingAction& action,
                                           Scene* scene, MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model),
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

void AtomMappingSceneTool::onLeftButtonPress(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonPress(event);
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(scene_pos, m_highlight_types);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        m_pressed_atom_item = atom_item;
    }
}

void AtomMappingSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonClick(event);
    if (m_action != MappingAction::REMOVE) {
        return;
    }

    if (m_pressed_atom_item == nullptr) {
        return;
    }
    auto* clicked_atom = m_pressed_atom_item->getAtom();
    int mapping_n = clicked_atom->getAtomMapNum();
    if (mapping_n < 1) {
        return;
    }
    // count the number of product atoms with the same mapping number as the
    // clicked atom
    unsigned int num_other_product_atoms = 0;
    std::unordered_set<const RDKit::Atom*> atoms_to_remove = {clicked_atom};
    auto mol = m_mol_model->getMol();
    for (auto* cur_atom : mol->atoms()) {
        if (cur_atom == clicked_atom ||
            cur_atom->getAtomMapNum() != mapping_n) {
            continue;
        }
        atoms_to_remove.insert(cur_atom);
        if (m_mol_model->isProductAtom(cur_atom)) {
            ++num_other_product_atoms;
        }
    }
    if (num_other_product_atoms > 0 &&
        m_mol_model->isProductAtom(clicked_atom)) {
        // there will still be remaining product atoms with this mapping
        // number, so only remove the mapping from the clicked atom
        atoms_to_remove = {clicked_atom};
    }
    m_mol_model->setAtomMapping(atoms_to_remove, 0);
}

bool AtomMappingSceneTool::isValidMappingPair(
    const AtomItem* const pressed_atom_item,
    const AtomItem* const hovered_atom_item)
{
    if (hovered_atom_item == nullptr ||
        hovered_atom_item == pressed_atom_item) {
        return false;
    }
    auto* pressed_atom = pressed_atom_item->getAtom();
    auto* hovered_atom = hovered_atom_item->getAtom();
    return (m_mol_model->isReactantAtom(pressed_atom) !=
            m_mol_model->isReactantAtom(hovered_atom));
}

void AtomMappingSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragMove(event);
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

    if (isValidMappingPair(m_pressed_atom_item, atom_item)) {
        m_release_atom_item = atom_item;
        // snap the line to the hovered atom
        scene_pos = atom_item->scenePos();
    } else {
        // the hovered atom item isn't valid for mapping, so ignore it
        atom_item = nullptr;
        m_release_atom_item = nullptr;
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

void AtomMappingSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragRelease(event);
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
    auto items = StandardSceneToolBase::getGraphicsItems();
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

QPixmap AtomMappingSceneTool::createDefaultCursorPixmap() const
{
    QString path = m_action == MappingAction::ADD
                       ? ":/icons/reaction_map_atoms.svg"
                       : ":/icons/reaction_unmap_atoms.svg";
    return cursor_hint_from_svg(path);
}

} // namespace sketcher
} // namespace schrodinger
