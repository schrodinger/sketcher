#include "schrodinger/sketcher/molviewer/scene_utils.h"

#include <QBitmap>
#include <QColor>
#include <QGraphicsItem>
#include <QPixmap>
#include <QPainter>
#include <QPen>
#include <QRect>
#include <QRegion>
#include <QSvgRenderer>

#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

std::tuple<std::vector<QGraphicsItem*>,
           std::unordered_map<const RDKit::Atom*, AtomItem*>,
           std::unordered_map<const RDKit::Bond*, BondItem*>>
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
        return {{}, {}, {}};
    }
    RDKit::Conformer conformer = shorten_attachment_point_bonds(mol);

    std::vector<QGraphicsItem*> all_items;
    std::unordered_map<const RDKit::Atom*, AtomItem*> atom_to_atom_item;
    std::unordered_map<const RDKit::Bond*, BondItem*> bond_to_bond_item;

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

    return {all_items, atom_to_atom_item, bond_to_bond_item};
}

void update_conf_for_mol_graphics_items(const QList<QGraphicsItem*>& atom_items,
                                        const QList<QGraphicsItem*>& bond_items,
                                        const RDKit::ROMol& mol)
{
    auto conf = mol.getConformer();
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

} // namespace sketcher
} // namespace schrodinger
