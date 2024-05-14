#include "schrodinger/sketcher/molviewer/scene_utils.h"

#include <QBitmap>
#include <QColor>
#include <QGraphicsItem>
#include <QLineF>
#include <QPixmap>
#include <QPainter>
#include <QPainterPath>
#include <QPen>
#include <QPointF>
#include <QRect>
#include <QRegion>
#include <QSvgRenderer>
#include <QTransform>

#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/sgroup.h"
#include "schrodinger/rdkit_extensions/variable_attachment_bond.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/sgroup_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

static QPainterPath
get_predictive_highlighting_path_for_bond(const RDKit::Bond* bond,
                                          const RDKit::Conformer& conf);

std::tuple<std::vector<QGraphicsItem*>,
           std::unordered_map<const RDKit::Atom*, AtomItem*>,
           std::unordered_map<const RDKit::Bond*, BondItem*>,
           std::unordered_map<const RDKit::SubstanceGroup*, SGroupItem*>>
create_graphics_items_for_mol(const RDKit::ROMol* mol, const Fonts& fonts,
                              AtomItemSettings& atom_item_settings,
                              BondItemSettings& bond_item_settings,
                              const bool draw_attachment_points)
{
    unsigned int num_atoms = mol->getNumAtoms();
    if (num_atoms == 0) {
        // If there are no atoms, then there's nothing more to do.  Also, if
        // there are no atoms, then shorten_attachment_point_bonds() will raise
        // a ConformerException.
        return {{}, {}, {}, {}};
    }
    RDKit::Conformer conformer =
        get_conformer_with_shortened_attachment_point_bonds(*mol);

    std::vector<QGraphicsItem*> all_items;
    std::unordered_map<const RDKit::Atom*, AtomItem*> atom_to_atom_item;
    std::unordered_map<const RDKit::Bond*, BondItem*> bond_to_bond_item;
    std::unordered_map<const RDKit::SubstanceGroup*, SGroupItem*>
        s_group_to_s_group_items;

    // create atom items
    for (std::size_t i = 0; i < num_atoms; ++i) {
        const auto* atom = mol->getAtomWithIdx(i);
        if (!draw_attachment_points && is_attachment_point(atom)) {
            continue;
        }
        const auto pos = conformer.getAtomPos(i);
        auto* atom_item = new AtomItem(atom, fonts, atom_item_settings);
        atom_item->setPos(to_scene_xy(pos));
        atom_to_atom_item[atom] = atom_item;
        if (rdkit_extensions::is_dummy_atom_for_variable_attachment_bond(
                atom)) {
            // hide the dummy atoms for variable attachment bonds
            atom_item->setVisible(false);
        }
        all_items.push_back(atom_item);
    }

    // create bond items
    for (auto bond : mol->bonds()) {
        if (!draw_attachment_points && is_attachment_point_bond(bond)) {
            continue;
        }
        const auto* from_atom_item = atom_to_atom_item[bond->getBeginAtom()];
        const auto* to_atom_item = atom_to_atom_item[bond->getEndAtom()];
        auto* bond_item = new BondItem(bond, *from_atom_item, *to_atom_item,
                                       fonts, bond_item_settings);
        bond_to_bond_item[bond] = bond_item;
        all_items.push_back(bond_item);
    }

    // create substance group items
    for (auto& sgroup : getSubstanceGroups(*mol)) {
        SGroupItem* sgroup_item = new SGroupItem(sgroup, fonts);
        s_group_to_s_group_items[&sgroup] = sgroup_item;
        all_items.push_back(sgroup_item);
    }

    return {all_items, atom_to_atom_item, bond_to_bond_item,
            s_group_to_s_group_items};
}

void update_conf_for_mol_graphics_items(
    const QList<QGraphicsItem*>& atom_items,
    const QList<QGraphicsItem*>& bond_items,
    const QList<QGraphicsItem*>& sgroup_items, const RDKit::ROMol& mol)
{
    auto conf = get_conformer_with_shortened_attachment_point_bonds(mol);
    for (auto* item : atom_items) {
        auto* atom_item = qgraphicsitem_cast<AtomItem*>(item);
        auto pos = conf.getAtomPos(atom_item->getAtom()->getIdx());
        atom_item->setPos(to_scene_xy(pos));
        atom_item->updateCachedData();
    }
    for (auto* item : bond_items) {
        auto* bond_item = qgraphicsitem_cast<BondItem*>(item);
        bond_item->updateCachedData();
    }
    for (auto* item : sgroup_items) {
        auto* sgroup_item = qgraphicsitem_cast<SGroupItem*>(item);
        sgroup_item->updateCachedData();
    }
}

QPixmap render_text_to_pixmap(const QString& text, const QFont& font,
                              const QColor& color)
{
    auto fm = QFontMetrics(font);
    auto bounding_rect = fm.boundingRect(text);
    // Add an extra pixel of width and height to the bounding rect so that we
    // don't cut off any of the text
    bounding_rect.adjust(0, 0, 1, 1);
    auto pixmap = QPixmap(bounding_rect.size());
    pixmap.fill(Qt::transparent);
    // pixmaps always start at (0, 0), so make sure the bounding rect starts
    // there too
    bounding_rect.moveTopLeft({0, 0});
    {
        auto painter = QPainter(&pixmap);
        painter.setFont(font);
        painter.setPen({color});
        painter.drawText(bounding_rect, Qt::AlignLeft | Qt::AlignTop, text);
    } // end the painter
    return pixmap;
}

QPixmap cursor_hint_from_svg(const QString& path, const bool recolor)
{
    QSvgRenderer renderer(path);
    renderer.setAspectRatioMode(Qt::KeepAspectRatio);
    auto size = renderer.defaultSize();
    size.scale(CURSOR_HINT_IMAGE_SIZE, CURSOR_HINT_IMAGE_SIZE,
               Qt::KeepAspectRatio);
    QPixmap pixmap(size);
    pixmap.fill(Qt::transparent);
    {
        // paint the SVG image to our pixmap
        QPainter painter(&pixmap);
        painter.setRenderHints(QPainter::Antialiasing |
                               QPainter::SmoothPixmapTransform);
        renderer.render(&painter, pixmap.rect());
    } // end the painter
    if (recolor) {
        // Create a mask denoting the dark gray pixels.  createMaskFromColor
        // requires that colors match exactly, so we first convert the image to
        // 8-bit color depth so that the mask won't miss pixels that are
        // imperceptibly different shades of gray.
        auto image = pixmap.toImage();
        image.convertTo(QImage::Format_Indexed8);
        auto image_mask =
            image.createMaskFromColor(TOOL_BUTTON_DARK_GRAY, Qt::MaskOutColor);
        auto mask = QBitmap::fromImage(image_mask);

        // Use the mask to paint CURSOR_HINT_COLOR over all gray pixels, but
        // keep the existing alpha channel data
        QPainter painter(&pixmap);
        painter.setPen(QColor(CURSOR_HINT_COLOR));
        painter.setRenderHints(QPainter::Antialiasing |
                               QPainter::SmoothPixmapTransform);
        painter.setCompositionMode(QPainter::CompositionMode_SourceIn);
        painter.drawPixmap(pixmap.rect(), mask);
    } // end the painter

    // Crop off the empty border of the image.  Otherwise the hint will wind up
    // too far away from the cursor
    auto bounding_rect = QRegion(pixmap.mask()).boundingRect();
    return pixmap.copy(bounding_rect);
}

QPixmap get_arrow_cursor_pixmap()
{
    QSvgRenderer renderer(ARROW_CURSOR_PATH);
    QPixmap pixmap(renderer.defaultSize());
    pixmap.fill(Qt::transparent);
    {
        // paint the SVG image to our pixmap
        QPainter painter(&pixmap);
        painter.setRenderHints(QPainter::Antialiasing |
                               QPainter::SmoothPixmapTransform);
        renderer.render(&painter, pixmap.rect());
    } // end the painter
    return pixmap;
}

QPainterPath get_predictive_highlighting_path_for_s_group_atoms_and_bonds(
    const RDKit::SubstanceGroup& s_group)
{
    QPainterPath path;
    path.setFillRule(Qt::WindingFill);
    const auto& mol = s_group.getOwningMol();
    const auto conf = get_conformer_with_shortened_attachment_point_bonds(mol);
    for (auto atom_idx : s_group.getAtoms()) {
        auto* atom = mol.getAtomWithIdx(atom_idx);
        auto cur_path = get_predictive_highlighting_path_for_atom(atom);
        auto& atom_pos = conf.getAtomPos(atom_idx);
        // translate the path out of the atom's local coordinate system
        cur_path.translate(to_scene_xy(atom_pos));
        path.addPath(cur_path);
    }
    for (auto* bond : rdkit_extensions::get_bonds_within_sgroup(s_group)) {
        auto cur_path = get_predictive_highlighting_path_for_bond(bond, conf);
        auto& bond_pos = conf.getAtomPos(bond->getBeginAtomIdx());
        // translate the path out of the bond's local coordinate system
        cur_path.translate(to_scene_xy(bond_pos));
        path.addPath(cur_path);
    }
    return path;
}

static QPainterPath get_highlighting_path_for_atom(const RDKit::Atom* atom,
                                                   const qreal radius)
{
    QPainterPath path;
    if (!is_attachment_point(atom)) {
        path.addEllipse(QPointF(0, 0), radius, radius);
        return path;
    } else {
        path.setFillRule(Qt::WindingFill);
        qreal squiggle_width = ATTACHMENT_POINT_SQUIGGLE_NUMBER_OF_WAVES *
                               ATTACHMENT_POINT_SQUIGGLE_WIDTH_PER_WAVE;
        qreal half_squiggle_width = 0.5 * squiggle_width;
        path.addRect(-half_squiggle_width, -radius, squiggle_width, 2 * radius);
        path.addEllipse(QPointF(-half_squiggle_width, 0.0), radius, radius);
        path.addEllipse(QPointF(half_squiggle_width, 0.0), radius, radius);

        qreal angle = get_attachment_point_line_angle(atom);
        QTransform transform;
        transform.rotate(angle);
        return transform.map(path);
    }
}

QPainterPath get_selection_highlighting_path_for_atom(const RDKit::Atom* atom)
{
    return get_highlighting_path_for_atom(atom,
                                          ATOM_SELECTION_HIGHLIGHTING_RADIUS);
}

QPainterPath get_predictive_highlighting_path_for_atom(const RDKit::Atom* atom)
{
    return get_highlighting_path_for_atom(atom,
                                          ATOM_PREDICTIVE_HIGHLIGHTING_RADIUS);
}

static QPainterPath get_highlighting_path_for_bond(const RDKit::Bond* bond,
                                                   const qreal half_width,
                                                   const RDKit::Conformer& conf)
{
    const auto& begin_pos = conf.getAtomPos(bond->getBeginAtomIdx());
    const auto& end_pos = conf.getAtomPos(bond->getEndAtomIdx());
    auto bond_vector = end_pos - begin_pos;
    QLineF bond_line = QLineF(QPointF(0, 0), to_scene_xy(bond_vector));
    return path_around_line(bond_line, half_width);
}

QPainterPath get_selection_highlighting_path_for_bond(const RDKit::Bond* bond)
{
    const auto conf = get_conformer_with_shortened_attachment_point_bonds(
        bond->getOwningMol());
    return get_highlighting_path_for_bond(
        bond, BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH, conf);
}

QPainterPath get_predictive_highlighting_path_for_bond(const RDKit::Bond* bond)
{
    const auto conf = get_conformer_with_shortened_attachment_point_bonds(
        bond->getOwningMol());
    return get_predictive_highlighting_path_for_bond(bond, conf);
}

/**
 * @overload This overload accepts a conformer where the attachment point bonds
 * have already been shortened.  This avoids the need to recalculate the
 * conformer on a per-bond basis, as that could lead to O(N^2) behavior.
 */
static QPainterPath
get_predictive_highlighting_path_for_bond(const RDKit::Bond* bond,
                                          const RDKit::Conformer& conf)
{
    auto path = get_highlighting_path_for_bond(
        bond, BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH, conf);
    auto atoms = rdkit_extensions::get_variable_attachment_atoms(bond);
    if (atoms.empty()) {
        // this isn't a variable attachment bond, so we don't need to add
        // anything else to the path
        return path;
    }

    // this is a variable attachment bond, so we want to add highlighting for
    // all of the variable attachment atoms and any bonds between those atoms
    path.setFillRule(Qt::WindingFill);
    std::unordered_set<RDKit::Bond*> added_bonds;
    const auto& mol = bond->getOwningMol();
    // we'll need to translate all of the paths below to use the same local
    // coordinate system as bond, so we fetch the bond's coordinates
    auto local_origin = conf.getAtomPos(bond->getBeginAtomIdx());
    for (auto* cur_atom : atoms) {
        auto atom_path = get_predictive_highlighting_path_for_atom(cur_atom);
        auto& atom_pos = conf.getAtomPos(cur_atom->getIdx());
        auto translate_by = to_scene_xy(atom_pos - local_origin);
        atom_path.translate(translate_by);
        path.addPath(atom_path);
        for (auto* cur_bond : mol.atomBonds(cur_atom)) {
            if (!added_bonds.count(cur_bond) &&
                atoms.count(cur_bond->getOtherAtom(cur_atom))) {
                // we haven't added a path for this bond yet and it's between
                // two variable attachment atoms, so we want to highlight it
                added_bonds.insert(cur_bond);
                auto bond_path = get_highlighting_path_for_bond(
                    cur_bond, BOND_PREDICTIVE_HIGHLIGHTING_HALF_WIDTH, conf);
                const auto& bond_pos =
                    conf.getAtomPos(cur_bond->getBeginAtomIdx());
                auto translate_by = to_scene_xy(bond_pos - local_origin);
                bond_path.translate(translate_by);
                path.addPath(bond_path);
            }
        }
    }
    return path;
}

QPainterPath path_around_line(const QLineF& line, const qreal half_width)
{
    QLineF normal = line.normalVector();
    normal.setLength(half_width);
    QPointF offset = normal.p2() - normal.p1();
    QPointF p1 = line.p1();
    QPointF p2 = line.p2();
    QPainterPath path;
    path.moveTo(p1 + offset);
    path.lineTo(p2 + offset);
    path.lineTo(p2 - offset);
    path.lineTo(p1 - offset);
    path.closeSubpath();
    return path;
}

} // namespace sketcher
} // namespace schrodinger
