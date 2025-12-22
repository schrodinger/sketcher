#include "schrodinger/sketcher/molviewer/scene_utils.h"

#include <QBitmap>
#include <QColor>
#include <QFile>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QLineF>
#include <QPixmap>
#include <QPainter>
#include <QPainterPath>
#include <QPen>
#include <QPointF>
#include <QRect>
#include <QRegularExpression>
#include <QRegion>
#include <QSvgRenderer>
#include <QTransform>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/rdkit/sgroup.h"
#include "schrodinger/sketcher/rdkit/variable_attachment_bond_core.h"
#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/amino_acid_item.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/chem_monomer_item.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/monomer_utils.h"
#include "schrodinger/sketcher/molviewer/non_molecular_item.h"
#include "schrodinger/sketcher/molviewer/nucleic_acid_base_item.h"
#include "schrodinger/sketcher/molviewer/nucleic_acid_phosphate_item.h"
#include "schrodinger/sketcher/molviewer/nucleic_acid_sugar_item.h"
#include "schrodinger/sketcher/molviewer/sgroup_item.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/rdkit/monomer_connectors.h"

namespace
{
const std::string PEPTIDE_POLYMER_PREFIX = "PEPTIDE";
// According to HELM, DNA is a subtype of RNA, so DNA also uses the RNA prefix
const std::string NUCLEOTIDE_POLYMER_PREFIX = "RNA";
} // namespace

namespace schrodinger
{
namespace sketcher
{

static QPainterPath
get_predictive_highlighting_path_for_bond(const RDKit::Bond* bond,
                                          const RDKit::Conformer& conf);

/**
 * Determine what type of monomer the given atom represents.
 */
MonomerType get_monomer_type(const RDKit::Atom* atom)
{
    const auto* monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom has no monomer info");
    }
    const auto* res_info =
        dynamic_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info);
    if (res_info == nullptr) {
        return MonomerType::CHEM;
    }
    const auto& chain_id = res_info->getChainId();
    if (chain_id.starts_with(PEPTIDE_POLYMER_PREFIX)) {
        return MonomerType::PEPTIDE;
    } else if (chain_id.starts_with(NUCLEOTIDE_POLYMER_PREFIX)) {
        const auto& res_name = res_info->getResidueName();
        if (res_name.empty()) {
            return MonomerType::NA_BASE;
        }
        auto last_char = std::tolower(res_name.back());
        if (last_char == 'p') {
            return MonomerType::NA_PHOSPHATE;
        } else if (last_char == 'r') {
            return MonomerType::NA_SUGAR;
        }
        return MonomerType::NA_BASE;
    }
    return MonomerType::CHEM;
}

/**
 * Construct and return a monomer graphics item for representing the given atom.
 */
AbstractMonomerItem* get_monomer_graphics_item(const RDKit::Atom* atom,
                                               const Fonts& fonts)
{
    switch (get_monomer_type(atom)) {
        case MonomerType::PEPTIDE:
            return new AminoAcidItem(atom, fonts);
        case MonomerType::CHEM:
            return new ChemMonomerItem(atom, fonts);
        case MonomerType::NA_BASE:
            return new NucleicAcidBaseItem(atom, fonts);
        case MonomerType::NA_PHOSPHATE:
            return new NucleicAcidPhosphateItem(atom, fonts);
        case MonomerType::NA_SUGAR:
            return new NucleicAcidSugarItem(atom, fonts);
        default:
            throw std::runtime_error("Unrecognized monomer type");
    }
}

std::tuple<std::vector<QGraphicsItem*>,
           std::unordered_map<const RDKit::Atom*, QGraphicsItem*>,
           std::unordered_map<const RDKit::Bond*, QGraphicsItem*>,
           std::unordered_map<const RDKit::Bond*, QGraphicsItem*>,
           std::unordered_map<const RDKit::SubstanceGroup*, SGroupItem*>>
create_graphics_items_for_mol(const RDKit::ROMol* mol, const Fonts& fonts,
                              const AtomDisplaySettings& atom_display_settings,
                              const BondDisplaySettings& bond_display_settings,
                              const bool draw_attachment_points)
{
    unsigned int num_atoms = mol->getNumAtoms();
    if (num_atoms == 0) {
        // If there are no atoms, then there's nothing more to do.  Also, if
        // there are no atoms, then shorten_attachment_point_bonds() will raise
        // a ConformerException.
        return {{}, {}, {}, {}, {}};
    }
    RDKit::Conformer conformer = mol->getConformer();

    std::vector<QGraphicsItem*> all_items;
    std::unordered_map<const RDKit::Atom*, QGraphicsItem*> atom_to_atom_item;
    std::unordered_map<const RDKit::Bond*, QGraphicsItem*> bond_to_bond_item;
    std::unordered_map<const RDKit::Bond*, QGraphicsItem*>
        bond_to_secondary_connection_item;
    std::unordered_map<const RDKit::SubstanceGroup*, SGroupItem*>
        s_group_to_s_group_items;

    // create atom items
    for (std::size_t i = 0; i < num_atoms; ++i) {
        const auto* atom = mol->getAtomWithIdx(i);
        if (!draw_attachment_points && is_attachment_point(atom)) {
            continue;
        }
        const auto pos = conformer.getAtomPos(i);
        QGraphicsItem* atom_item =
            is_atom_monomeric(atom)
                ? static_cast<QGraphicsItem*>(
                      get_monomer_graphics_item(atom, fonts))
                : new AtomItem(atom, fonts, atom_display_settings);
        atom_item->setPos(to_scene_xy(pos));
        atom_to_atom_item[atom] = atom_item;
        if (is_dummy_atom_for_variable_attachment_bond(atom)) {
            // hide the dummy atoms for variable attachment bonds
            atom_item->setVisible(false);
        }
        all_items.push_back(atom_item);
    }

    // create bond items
    for (auto* bond : mol->bonds()) {
        if (!draw_attachment_points && is_attachment_point_bond(bond)) {
            continue;
        }
        const auto* from_atom = bond->getBeginAtom();
        const auto* to_atom = bond->getEndAtom();
        const auto* from_graphics_item = atom_to_atom_item[from_atom];
        const auto* to_graphics_item = atom_to_atom_item[to_atom];
        const auto* from_atom_item =
            qgraphicsitem_cast<const AtomItem*>(from_graphics_item);
        const auto* to_atom_item =
            qgraphicsitem_cast<const AtomItem*>(to_graphics_item);
        const auto* from_monomer_item =
            dynamic_cast<const AbstractMonomerItem*>(from_graphics_item);
        const auto* to_monomer_item =
            dynamic_cast<const AbstractMonomerItem*>(to_graphics_item);
        if (from_atom_item != nullptr && to_atom_item != nullptr) {
            auto bond_item = new BondItem(bond, *from_atom_item, *to_atom_item,
                                          fonts, bond_display_settings);
            bond_to_bond_item[bond] = bond_item;
            all_items.push_back(bond_item);
        } else if (from_monomer_item != nullptr && to_monomer_item != nullptr) {
            auto connector_item = new MonomerConnectorItem(
                bond, *from_monomer_item, *to_monomer_item);
            bond_to_bond_item[bond] = connector_item;
            all_items.push_back(connector_item);
            if (contains_two_monomer_linkages(bond)) {
                auto secondary_connector_item = new MonomerConnectorItem(
                    bond, *from_monomer_item, *to_monomer_item,
                    /* is_secondary_connection = */ true);
                bond_to_secondary_connection_item[bond] =
                    secondary_connector_item;
                all_items.push_back(secondary_connector_item);
            }
        }
        // we skip bonds that go between a monomer and an atomistic atom
    }

    // create substance group items
    for (auto& sgroup : getSubstanceGroups(*mol)) {
        if (rdkit_extensions::is_polymer_annotation_s_group(sgroup) ||
            rdkit_extensions::is_supplementary_information_s_group(sgroup)) {
            // this isn't an actual S-group; it's just additional data about a
            // monomeric model
            continue;
        }
        SGroupItem* sgroup_item = new SGroupItem(sgroup, fonts);
        s_group_to_s_group_items[&sgroup] = sgroup_item;
        all_items.push_back(sgroup_item);
    }

    return {all_items, atom_to_atom_item, bond_to_bond_item,
            bond_to_secondary_connection_item, s_group_to_s_group_items};
}

void update_conf_for_mol_graphics_items(
    const QList<QGraphicsItem*>& atom_items,
    const QList<QGraphicsItem*>& bond_items,
    const QList<QGraphicsItem*>& sgroup_items, const RDKit::ROMol& mol)
{
    auto conf = mol.getConformer();
    for (auto* item : atom_items) {
        auto* atom_item = static_cast<AbstractAtomOrMonomerItem*>(item);
        auto pos = conf.getAtomPos(atom_item->getAtom()->getIdx());
        atom_item->setPos(to_scene_xy(pos));
        atom_item->updateCachedData();
    }
    for (auto* item : bond_items) {
        auto* bond_item = static_cast<AbstractBondOrConnectorItem*>(item);
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

QPixmap cursor_hint_from_graphics_item(QGraphicsItem* graphics_item,
                                       const qreal min_scene_size)
{
    QGraphicsScene scene;
    QPixmap pixmap(CURSOR_HINT_IMAGE_SIZE, CURSOR_HINT_IMAGE_SIZE);
    scene.addItem(graphics_item);
    QRectF render_source = scene.sceneRect();
    render_source.setWidth(std::max(render_source.width(), min_scene_size));
    render_source.setHeight(std::max(render_source.height(), min_scene_size));

    pixmap.fill(Qt::transparent);
    {
        QPainter painter(&pixmap);
        painter.setRenderHints(QPainter::Antialiasing |
                               QPainter::SmoothPixmapTransform);
        scene.render(&painter, QRectF(), render_source);
    }

    // remove graphics_item from the scene, otherwise it'll get destroyed by Qt
    // when the scene is destroyed at the end of this function
    scene.removeItem(graphics_item);
    return pixmap;
}

QPixmap get_arrow_cursor_pixmap(const QColor& arrow_color,
                                const QColor& outline_color)
{
    // read in the contents of the SVG
    QString svg_contents;
    {
        QFile file(ARROW_CURSOR_PATH);
        bool success = file.open(QIODevice::ReadOnly);
        if (!success) {
            throw std::runtime_error("Could not open cursor pixmap file");
        }
        svg_contents = file.readAll();
    } // close the file

    // replace the colors in the SVG with the colors that were passed in to this
    // function
    QRegularExpression fill_color_re(R"(fill\s*:\s*#\w{3,8})");
    svg_contents.replace(fill_color_re, QString("fill:") + arrow_color.name());
    QRegularExpression stroke_color_re(R"(stroke\s*:\s*#\w{3,8})");
    svg_contents.replace(stroke_color_re,
                         QString("stroke:") + outline_color.name());

    QSvgRenderer renderer(svg_contents.toUtf8());
    // figure out how large the arrow should be, taking CURSOR_SCALE into
    // account
    auto scaled_box = renderer.viewBoxF();
    scaled_box.setWidth(scaled_box.width() * CURSOR_SCALE);
    scaled_box.setHeight(scaled_box.height() * CURSOR_SCALE);
    // shift the image so that the tip of the arrow lines up exactly with the
    // hotspot coordinates.  We can't scale the hotspot coordinates themselves
    // (at least not precisely) because the coordinates must be integers.
    scaled_box.translate(CURSOR_HOTSPOT_X - CURSOR_SCALE,
                         CURSOR_HOTSPOT_Y - CURSOR_SCALE);

    // render the svg
    QPixmap pixmap(scaled_box.toAlignedRect().size());
    pixmap.fill(Qt::transparent);
    {
        // paint the SVG image to our pixmap
        QPainter painter(&pixmap);
        painter.setRenderHints(QPainter::Antialiasing |
                               QPainter::SmoothPixmapTransform);
        renderer.render(&painter, scaled_box);
    } // end the painter
    return pixmap;
}

QPainterPath get_predictive_highlighting_path_for_s_group_atoms_and_bonds(
    const RDKit::SubstanceGroup& s_group)
{
    QPainterPath path;
    path.setFillRule(Qt::WindingFill);
    const auto& mol = s_group.getOwningMol();
    const auto conf = mol.getConformer();
    for (auto atom_idx : s_group.getAtoms()) {
        auto* atom = mol.getAtomWithIdx(atom_idx);
        auto cur_path = get_predictive_highlighting_path_for_atom(atom);
        auto& atom_pos = conf.getAtomPos(atom_idx);
        // translate the path out of the atom's local coordinate system
        cur_path.translate(to_scene_xy(atom_pos));
        path.addPath(cur_path);
    }
    for (auto* bond : get_bonds_within_sgroup(s_group)) {
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
    const auto conf = bond->getOwningMol().getConformer();
    return get_highlighting_path_for_bond(
        bond, BOND_SELECTION_HIGHLIGHTING_HALF_WIDTH, conf);
}

QPainterPath get_predictive_highlighting_path_for_bond(const RDKit::Bond* bond)
{
    const auto conf = bond->getOwningMol().getConformer();
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
    auto atoms = get_variable_attachment_atoms(bond);
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

bool item_matches_type_flag(QGraphicsItem* item,
                            InteractiveItemFlagType type_flag)
{
    // a mapping of graphics items types to InteractiveItemFlag for graphics
    // items that represent a single InteractiveItemFlag value
    static const std::unordered_map<int, InteractiveItemFlagType>
        type_to_interactive_item_flag{
            {SGroupItem::Type, InteractiveItemFlag::S_GROUP},
            {NonMolecularItem::Type, InteractiveItemFlag::NON_MOLECULAR},
            {AminoAcidItem::Type, InteractiveItemFlag::AA_MONOMER},
            {NucleicAcidPhosphateItem::Type, InteractiveItemFlag::NA_PHOSPHATE},
            {NucleicAcidSugarItem::Type, InteractiveItemFlag::NA_SUGAR},
            {NucleicAcidBaseItem::Type, InteractiveItemFlag::NA_BASE},
            {ChemMonomerItem::Type, InteractiveItemFlag::CHEM_MONOMER},
            {MonomerConnectorItem::Type,
             InteractiveItemFlag::MONOMER_CONNECTOR},
        };

    auto type = item->type();

    // AtomItems and BondItems can represent multiple InteractiveItemFlag values
    // (e.g. R-groups, attachment points, attachment point bonds), so we
    // consider those separately
    if (type == AtomItem::Type) {
        if ((type_flag & InteractiveItemFlag::ATOM) ==
            InteractiveItemFlag::ATOM) {
            // we want all atoms, regardless of whether or not they're
            // attachment points
            return true;
        } else if (type_flag & InteractiveItemFlag::ATOM) {
            auto* atom = static_cast<AtomItem*>(item)->getAtom();
            if (is_r_group(atom)) {
                return type_flag & InteractiveItemFlag::R_GROUP;
            } else if (is_attachment_point(atom)) {
                return type_flag & InteractiveItemFlag::ATTACHMENT_POINT;
            } else {
                return type_flag & InteractiveItemFlag::ATOM_NOT_R_NOT_AP;
            }
        }
    } else if (type == BondItem::Type) {
        if ((type_flag & InteractiveItemFlag::BOND) ==
            InteractiveItemFlag::BOND) {
            // we want all bonds, regardless of whether or not they're
            // attachment point bonds
            return true;
        } else if (type_flag & InteractiveItemFlag::BOND) {
            auto* bond = static_cast<BondItem*>(item)->getBond();
            return is_attachment_point_bond(bond) ==
                   static_cast<bool>(
                       type_flag & InteractiveItemFlag::ATTACHMENT_POINT_BOND);
        }
    }

    auto search = type_to_interactive_item_flag.find(type);
    if (search != type_to_interactive_item_flag.end()) {
        return type_flag & search->second;
    }
    return false;
}

template <typename T>
std::tuple<std::unordered_set<const RDKit::Atom*>,
           std::unordered_set<const RDKit::Bond*>,
           std::unordered_set<const RDKit::Bond*>,
           std::unordered_set<const RDKit::SubstanceGroup*>,
           std::unordered_set<const NonMolecularObject*>>
get_model_objects_for_graphics_items(const T& items)
{
    std::unordered_set<const RDKit::Atom*> atoms;
    std::unordered_set<const RDKit::Bond*> bonds;
    std::unordered_set<const RDKit::Bond*> secondary_connections;
    std::unordered_set<const RDKit::SubstanceGroup*> s_groups;
    std::unordered_set<const NonMolecularObject*> non_molecular_objects;
    for (auto item : items) {
        if (item_matches_type_flag(item,
                                   InteractiveItemFlag::ATOM_OR_MONOMER)) {
            const auto* atom_or_monomer_item =
                static_cast<const AbstractAtomOrMonomerItem*>(item);
            atoms.insert(atom_or_monomer_item->getAtom());
        } else if (item_matches_type_flag(item, InteractiveItemFlag::BOND)) {
            auto* bond_item = qgraphicsitem_cast<BondItem*>(item);
            bonds.insert(bond_item->getBond());
        } else if (item_matches_type_flag(
                       item, InteractiveItemFlag::MONOMER_CONNECTOR)) {
            auto* connector_item =
                qgraphicsitem_cast<MonomerConnectorItem*>(item);
            auto* bond = connector_item->getBond();
            if (connector_item->isSecondaryConnection()) {
                secondary_connections.insert(bond);
            } else {
                bonds.insert(bond);
            }
        } else if (auto s_group_item = qgraphicsitem_cast<SGroupItem*>(item)) {
            s_groups.insert(s_group_item->getSubstanceGroup());
        } else if (auto nonmolecular_item =
                       qgraphicsitem_cast<NonMolecularItem*>(item)) {
            non_molecular_objects.insert(
                nonmolecular_item->getNonMolecularObject());
        }
    }
    return {atoms, bonds, secondary_connections, s_groups,
            non_molecular_objects};
}

template SKETCHER_API
    std::tuple<std::unordered_set<const RDKit::Atom*>,
               std::unordered_set<const RDKit::Bond*>,
               std::unordered_set<const RDKit::Bond*>,
               std::unordered_set<const RDKit::SubstanceGroup*>,
               std::unordered_set<const NonMolecularObject*>>
    get_model_objects_for_graphics_items<QList<QGraphicsItem*>>(
        const QList<QGraphicsItem*>& items);

template SKETCHER_API
    std::tuple<std::unordered_set<const RDKit::Atom*>,
               std::unordered_set<const RDKit::Bond*>,
               std::unordered_set<const RDKit::Bond*>,
               std::unordered_set<const RDKit::SubstanceGroup*>,
               std::unordered_set<const NonMolecularObject*>>
    get_model_objects_for_graphics_items<std::unordered_set<QGraphicsItem*>>(
        const std::unordered_set<QGraphicsItem*>& items);

} // namespace sketcher
} // namespace schrodinger
