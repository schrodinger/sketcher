#include "schrodinger/sketcher/tool/attachment_point_scene_tool.h"

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

HintSquiggleItem::HintSquiggleItem(QGraphicsItem* parent) :
    QGraphicsPathItem(parent)
{
    setZValue(static_cast<qreal>(ZOrder::HINT));
    setPen(QPen(STRUCTURE_HINT_COLOR));
    setVisible(false);
}

DrawAttachmentPointSceneTool::DrawAttachmentPointSceneTool(
    Scene* scene, MolModel* mol_model) :
    AbstractDrawSceneTool(scene, mol_model)
{
}

std::vector<QGraphicsItem*> DrawAttachmentPointSceneTool::getGraphicsItems()
{
    auto items = AbstractDrawSceneTool::getGraphicsItems();
    items.push_back(&m_hint_squiggle_item);
    return items;
}

QPixmap DrawAttachmentPointSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/enumeration_attachment_point.svg");
}

void DrawAttachmentPointSceneTool::onBondClicked(const RDKit::Bond* const bond)
{
    cycleBond(bond);
}

bool DrawAttachmentPointSceneTool::shouldDrawBondForClickOnAtom(
    const RDKit::Atom* const atom) const
{
    return number_of_bound_attachment_points(atom) <
           MAX_ATTACHMENT_POINTS_PER_ATOM;
}

void DrawAttachmentPointSceneTool::addAtom(const RDGeom::Point3D& pos,
                                           const RDKit::Atom* const bound_to)
{
    auto& conf = m_mol_model->getMol()->getConformer();
    RDGeom::Point3D bound_to_pos = conf.getAtomPos(bound_to->getIdx());
    RDGeom::Point3D bond_offset = (pos - bound_to_pos);
    RDGeom::Point3D new_pos = bound_to_pos + bond_offset;
    m_mol_model->addAttachmentPoint(new_pos, bound_to);
}

std::pair<QPointF, const RDKit::Atom*>
DrawAttachmentPointSceneTool::getBondEndInMousedDirection(
    const QPointF& start, const RDKit::Atom* const start_atom,
    const QPointF& mouse_pos) const
{
    QPointF bond_offset =
        getDefaultBondOffsetInMousedDirection(start, mouse_pos);
    QPointF bond_end = start + bond_offset;
    return {bond_end, nullptr};
}

std::tuple<RDGeom::Point3D, QPointF, const RDKit::Atom*>
DrawAttachmentPointSceneTool::getDefaultBondPosition(
    const RDKit::Atom* const atom) const
{
    auto [bond_end, scene_bond_end] = getInitialDefaultBondPosition(atom);
    return {bond_end, scene_bond_end, nullptr};
}

std::tuple<bool, QPointF, const RDKit::Atom*>
DrawAttachmentPointSceneTool::getDragStartInfo() const
{
    auto* item = m_scene->getTopInteractiveItemAt(
        m_mouse_press_scene_pos, InteractiveItemFlag::MOLECULAR_NOT_AP);
    if (auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        const auto* atom = atom_item->getAtom();
        if (number_of_bound_attachment_points(atom) <
            MAX_ATTACHMENT_POINTS_PER_ATOM) {
            return {true, atom_item->pos(), atom};
        }
    }
    return {false, QPointF(), nullptr};
}

void DrawAttachmentPointSceneTool::setHintBondVisible(bool visible)
{
    AbstractDrawSceneTool::setHintBondVisible(visible);
    m_hint_squiggle_item.setVisible(visible);
}

void DrawAttachmentPointSceneTool::updateHintBondPath(const QPointF& start,
                                                      const QPointF& end)
{
    AbstractDrawSceneTool::updateHintBondPath(start, end);
    QLineF bond_line(start, end);
    qreal angle = -bond_line.normalVector().angle();
    auto path = get_wavy_line_path(ATTACHMENT_POINT_SQUIGGLE_NUMBER_OF_WAVES,
                                   ATTACHMENT_POINT_SQUIGGLE_WIDTH_PER_WAVE,
                                   ATTACHMENT_POINT_SQUIGGLE_HEIGHT, angle);
    QTransform transform;
    transform.translate(end.x(), end.y());
    auto translated_path = transform.map(path);
    m_hint_squiggle_item.setPath(translated_path);
}

} // namespace sketcher
} // namespace schrodinger
