#include "schrodinger/sketcher/tool/draw_chain_scene_tool.h"

#include <algorithm>

#include <QFont>
#include <QLineF>
#include <QPolygonF>

#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"

namespace schrodinger
{
namespace sketcher
{

HintChainItem::HintChainItem(QGraphicsItem* parent) : QGraphicsItemGroup(parent)
{
    setZValue(static_cast<qreal>(ZOrder::HINT));
    setVisible(false);

    m_bonds_item.setPen(QPen(STRUCTURE_HINT_COLOR));
    m_label_item.setFont(QFont(FONT_NAME, DRAW_CHAIN_FONT_SIZE));

    addToGroup(&m_bonds_item);
    addToGroup(&m_label_item);
}

void HintChainItem::setCoords(QList<QPointF> coords)
{
    QPainterPath path;
    auto polygon = QPolygonF(coords);
    path.addPolygon(polygon);
    m_bonds_item.setPath(path);

    m_label_item.setPos(coords.back());
    // The label displays the number of bonds, which is one less than the number
    // of atoms
    m_label_item.setText(QString::number(coords.size() - 1));
}

DrawChainSceneTool::DrawChainSceneTool(Scene* scene, MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model)
{
    m_highlight_types = InteractiveItemFlag::ATOM_NOT_AP;
}

std::vector<QGraphicsItem*> DrawChainSceneTool::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
    items.push_back(&m_hint_chain_item);
    return items;
}

void DrawChainSceneTool::onLeftButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragStart(event);
    m_hint_chain_item.show();
}

void DrawChainSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragMove(event);
    auto [start_pos, start_atom] = getStartPosAndAtom();
    auto atom_coords = get_bond_chain_atom_coords(start_pos, event->scenePos());
    m_hint_chain_item.setCoords(atom_coords);
}

void DrawChainSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragRelease(event);
    m_hint_chain_item.hide();

    auto [start_pos, start_atom] = getStartPosAndAtom();
    auto atom_coords_scene =
        get_bond_chain_atom_coords(start_pos, event->scenePos());
    std::vector<RDGeom::Point3D> atom_coords_mol;
    unsigned int skip_first_coord = start_atom == nullptr ? 0 : 1;
    std::transform(atom_coords_scene.begin() + skip_first_coord,
                   atom_coords_scene.end(), std::back_inserter(atom_coords_mol),
                   to_mol_xy);
    m_mol_model->addAtomChain(Element::C, atom_coords_mol, start_atom);
}

std::pair<QPointF, const RDKit::Atom*> DrawChainSceneTool::getStartPosAndAtom()
{
    QPointF start_pos;
    const RDKit::Atom* atom = nullptr;
    auto* item = m_scene->getTopInteractiveItemAt(
        m_mouse_press_scene_pos, InteractiveItemFlag::ATOM_NOT_AP);
    if (item == nullptr) {
        start_pos = m_mouse_press_scene_pos;
    } else {
        auto* atom_item = static_cast<AtomItem*>(item);
        start_pos = atom_item->pos();
        atom = atom_item->getAtom();
    }
    return {start_pos, atom};
}

QPixmap DrawChainSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(":/icons/bond_chain.svg");
}

QList<QPointF> get_bond_chain_atom_coords(QPointF start, QPointF end)
{
    // round the angle to the nearest pi/6 radians (30 degrees).  The comments
    // below will refer to the vector V, which is an infinite vector emerging
    // from the start point at this rounded angle.
    qreal angle = get_rounded_angle_radians(start, end);

    // Note that the first bond makes an angle of pi/6 radians (30 degrees) with
    // the vector V.  The second bond makes an angle of -pi/6 radians (-30
    // degrees), so the end point of the second bond is *on* vector V.

    // Calculate the length of a bond projected onto the vector V.
    const qreal LENGTH_OF_PROJECTED_BOND = BOND_LENGTH * qCos(M_PI / 6.0);

    // create a QPointF that represents the end point of the bond projected onto
    // vector V
    auto projected_bond_onto_vector_mol =
        RDGeom::Point3D(qCos(angle), qSin(angle), 0);
    projected_bond_onto_vector_mol *= LENGTH_OF_PROJECTED_BOND;
    auto projected_bond_onto_vector_scene =
        to_scene_xy(projected_bond_onto_vector_mol);

    // Create a QPointF that represents the perpendicular from the vector V to
    // the end point of the first bond.  (In other words, start +
    // projected_bond_onto_vector_scene + perpendicular_scene gives the end
    // point of the first bond.)
    auto projected_bond_line =
        QLineF(QPointF(0, 0), projected_bond_onto_vector_scene);
    auto perpendicular_scene =
        projected_bond_line.unitVector().normalVector().p2();
    constexpr qreal SIN_PI_OVER_6 = 0.5;
    perpendicular_scene *= VIEW_SCALE * BOND_LENGTH * SIN_PI_OVER_6;

    // figure out how many bonds we need
    qreal length = QLineF(start, end).length();
    int num_bonds =
        std::round(length / (VIEW_SCALE * LENGTH_OF_PROJECTED_BOND));
    if (num_bonds == 0) {
        // Always include at least one bond
        num_bonds = 1;
    }

    // actually generate the coordinates
    QList<QPointF> coords;
    for (int i = 0; i <= num_bonds; i++) {
        auto cur_coord = start + (i * projected_bond_onto_vector_scene) +
                         ((i % 2) * perpendicular_scene);
        coords.push_back(cur_coord);
    }
    return coords;
}

} // namespace sketcher
} // namespace schrodinger
