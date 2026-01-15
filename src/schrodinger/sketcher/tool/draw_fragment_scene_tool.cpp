#include "schrodinger/sketcher/tool/draw_fragment_scene_tool.h"

#include <stdexcept>

#include <boost/shared_ptr.hpp>

#include <QGraphicsScene>
#include <QPainter>
#include <QPixmap>
#include <QPointF>

#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/sgroup_item.h"
#include "schrodinger/sketcher/rdkit/fragment.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

std::shared_ptr<DrawFragmentSceneTool>
get_draw_fragment_scene_tool(const RingTool& ring_tool, const Fonts& fonts,
                             const AtomDisplaySettings& atom_display_settings,
                             const BondDisplaySettings& bond_display_settings,
                             Scene* scene, MolModel* mol_model)
{
    std::string smiles = get_smiles_for_ring_tool(ring_tool);
    return get_draw_fragment_scene_tool(smiles, fonts, atom_display_settings,
                                        bond_display_settings, scene, mol_model,
                                        Format::EXTENDED_SMILES);
}

std::shared_ptr<DrawFragmentSceneTool>
get_draw_fragment_scene_tool(const std::string& text_mol, const Fonts& fonts,
                             const AtomDisplaySettings& atom_display_settings,
                             const BondDisplaySettings& bond_display_settings,
                             Scene* scene, MolModel* mol_model,
                             const Format format)
{
    auto mol = rdkit_extensions::to_rdkit(text_mol, format);
    prepare_mol(*mol);
    return std::make_shared<DrawFragmentSceneTool>(
        *mol, fonts, atom_display_settings, bond_display_settings, scene,
        mol_model);
}

std::string get_smiles_for_ring_tool(const RingTool& ring_tool)
{
    // CXSMILES coordinates ensure that the ring orientation matches the
    // orientation of the image on each toolbar button
    switch (ring_tool) {
        case (RingTool::CYCLOPROPANE):
            return "*C1CC1 "
                   "|(0.00439306,1.77451,;0.000679609,0.274518,;-0.752534,-1."
                   "02266,;0.747462,-1.02637,),$_AP1;;;$|";
        case (RingTool::CYCLOBUTANE):
            return "*C1CCC1 "
                   "|(-1.44716,1.4499,;-0.3875,0.388235,;1.1125,0.386814,;1."
                   "11108,-1.11318,;-0.388921,-1.11176,),$_AP1;;;;$|";
        case (RingTool::CYCLOPENTANE):
            return "*C1CCCC1 "
                   "|(-0.00118323,2.31331,;-0.000416,0.813313,;-1.21349,-0."
                   "0689851,;-0.749235,-1.49533,;0.750765,-1.49457,;1.21356,-0."
                   "0677437,),$_AP1;;;;;$|";
        case (RingTool::CYCLOPENTADIENE):
            return "*C1C=CC=C1 "
                   "|(-0.221326,3.16294,;-0.2265,1.66295,;0.983978,0.777089,;0."
                   "515536,-0.647889,;-0.984455,-0.642716,;-1.44306,0.785459,),"
                   "$_AP1;;;;;$|";
        case (RingTool::CYCLOHEXANE):
            return "*C1CCCCC1 "
                   "|(-0.0278384,2.57128,;-0.0115993,1.07137,;1.29548,0.335473,"
                   ";1.31172,-1.16444,;0.0208788,-1.92846,;-1.2862,-1.19257,;-"
                   "1.30244,0.307346,),$_AP1;;;;;;$|";
        case (RingTool::BENZENE):
            return "*C1=CC=CC=C1 "
                   "|(-0.0278384,2.57128,;-0.0115993,1.07137,;1.29548,0.335473,"
                   ";1.31172,-1.16444,;0.0208788,-1.92846,;-1.2862,-1.19257,;-"
                   "1.30244,0.307346,),$_AP1;;;;;;$|";
        case (RingTool::CYCLOHEPTANE):
            return "*C1CCCCCC1 "
                   "|(0.00929483,2.82499,;0.00435952,1.325,;-1.34923,0.67862,;-"
                   "1.68782,-0.782667,;-0.756448,-1.95849,;0.743544,-1.96342,;"
                   "1.68263,-0.793756,;1.35366,0.669727,),$_AP1;;;;;;;$|";
        case (RingTool::CYCLOOCTANE):
            return "*C1CCCCCCC1 "
                   "|(-1.17747,2.84108,;-0.603174,1.45538,;0.896826,1.45567,;1."
                   "9577,0.395223,;1.95799,-1.10478,;0.897542,-2.16565,;-0."
                   "602457,-2.16594,;-1.66333,-1.10549,;-1.66362,0.394506,),$_"
                   "AP1;;;;;;;;$|";
    }
    throw std::runtime_error("Unrecognized ring type");
}

AbstractHintFragmentItem::AbstractHintFragmentItem(
    const RDKit::ROMol& fragment, const Fonts& fonts,
    const AtomDisplaySettings& atom_display_settings,
    const BondDisplaySettings& bond_display_settings, QGraphicsItem* parent) :
    QGraphicsItemGroup(parent),
    m_frag(fragment),
    m_fonts(&fonts),
    m_atom_display_settings(atom_display_settings),
    m_bond_display_settings(bond_display_settings)
{
}

AbstractHintFragmentItem::AbstractHintFragmentItem(
    const AbstractHintFragmentItem& frag_item) :
    AbstractHintFragmentItem(frag_item.m_frag, *frag_item.m_fonts,
                             frag_item.m_atom_display_settings,
                             frag_item.m_bond_display_settings, nullptr)
{
}

void AbstractHintFragmentItem::updateSettings()
{
    m_atom_display_settings.setMonochromeColorScheme(STRUCTURE_HINT_COLOR);
    m_bond_display_settings.setMonochromeColorScheme(STRUCTURE_HINT_COLOR);
    // Disable atom indices for hints/previews
    m_atom_display_settings.m_show_atom_indices = false;
}

void AbstractHintFragmentItem::createGraphicsItems()
{
    auto [all_items, atom_to_atom_item, bond_to_bond_item,
          bond_to_secondary_connection_item, s_group_to_s_group_item] =
        create_graphics_items_for_mol(
            &m_frag, *m_fonts, m_atom_display_settings, m_bond_display_settings,
            /*draw_attachment_points = */ false);
    for (auto* item : all_items) {
        addToGroup(item);
    }
    for (auto& kv : atom_to_atom_item) {
        m_atom_items.append(kv.second);
    }
    for (auto& kv : bond_to_bond_item) {
        m_bond_items.append(kv.second);
    }
    for (auto& kv : bond_to_secondary_connection_item) {
        m_bond_items.append(kv.second);
    }
    for (auto& kv : s_group_to_s_group_item) {
        m_s_group_items.append(kv.second);
    }
}

HintFragmentItem::HintFragmentItem(
    const RDKit::ROMol& fragment, const Fonts& fonts,
    const AtomDisplaySettings& atom_display_settings,
    const BondDisplaySettings& bond_display_settings, QGraphicsItem* parent) :
    AbstractHintFragmentItem(fragment, fonts, atom_display_settings,
                             bond_display_settings, parent)
{
    std::transform(m_frag.bonds().begin(), m_frag.bonds().end(),
                   std::back_inserter(m_orig_bond_types),
                   [](const auto* bond) { return bond->getBondType(); });

    // We'll set this to visible once the mouse is over the scene
    setVisible(false);
    setZValue(static_cast<qreal>(ZOrder::HINT));

    updateSettings();
    createGraphicsItems();
}

void HintFragmentItem::updateSettings()
{
    AbstractHintFragmentItem::updateSettings();
    m_bond_display_settings.m_bond_width = FRAGMENT_HINT_BOND_WIDTH;
}

void HintFragmentItem::updateConformer(const RDKit::Conformer& conformer)
{
    m_frag.getConformer() = conformer;
    update_conf_for_mol_graphics_items(m_atom_items, m_bond_items,
                                       m_s_group_items, m_frag);
}

void HintFragmentItem::updateSingleBondMutations(
    const std::unordered_set<RDKit::Bond*>& bonds_to_mutate)
{
    // this class stores a copy of the molecule, so we need to use bond indices
    // rather than bond pointers
    std::vector<unsigned int> bond_idxs_to_mutate;
    std::transform(bonds_to_mutate.begin(), bonds_to_mutate.end(),
                   std::back_inserter(bond_idxs_to_mutate),
                   [](const auto* bond) { return bond->getIdx(); });
    std::sort(bond_idxs_to_mutate.begin(), bond_idxs_to_mutate.end());
    auto bonds_idxs_to_mutate_it = bond_idxs_to_mutate.begin();
    for (unsigned int bond_idx = 0; bond_idx < m_frag.getNumBonds();
         ++bond_idx) {
        RDKit::Bond::Bond::BondType cur_bond_type =
            RDKit::Bond::Bond::BondType::SINGLE;
        if (bonds_idxs_to_mutate_it != bond_idxs_to_mutate.end() &&
            bond_idx == *bonds_idxs_to_mutate_it) {
            ++bonds_idxs_to_mutate_it;
        } else {
            cur_bond_type = m_orig_bond_types[bond_idx];
        }
        m_frag.getBondWithIdx(bond_idx)->setBondType(cur_bond_type);
    }
    for (auto cur_bond_item : m_bond_items) {
        static_cast<BondItem*>(cur_bond_item)->updateCachedData();
    }
}

HintPixmapFragmentItem::HintPixmapFragmentItem(
    const RDKit::ROMol& fragment, const Fonts& fonts,
    const AtomDisplaySettings& atom_display_settings,
    const BondDisplaySettings& bond_display_settings, QGraphicsItem* parent) :
    AbstractHintFragmentItem(fragment, fonts, atom_display_settings,
                             bond_display_settings, parent)
{
    updateSettings();
    createGraphicsItems();
}

HintPixmapFragmentItem::HintPixmapFragmentItem(
    const AbstractHintFragmentItem& frag_item) :
    AbstractHintFragmentItem(frag_item)
{
    updateSettings();
    createGraphicsItems();
}

void HintPixmapFragmentItem::updateSettings()
{
    AbstractHintFragmentItem::updateSettings();
    m_bond_display_settings.m_bond_width = FRAGMENT_CURSOR_HINT_BOND_WIDTH;
    m_bond_display_settings.m_double_bond_spacing =
        FRAGMENT_CURSOR_HINT_BOND_SPACING;
}

DrawFragmentSceneTool::DrawFragmentSceneTool(
    const RDKit::ROMol& fragment, const Fonts& fonts,
    const AtomDisplaySettings& atom_display_settings,
    const BondDisplaySettings& bond_display_settings, Scene* scene,
    MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model),
    m_frag(fragment),
    m_hint_item(m_frag, fonts, atom_display_settings, bond_display_settings)
{
    m_highlight_types = InteractiveItemFlag::MOLECULAR_NOT_AP;
}

std::vector<QGraphicsItem*> DrawFragmentSceneTool::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
    items.push_back(&m_hint_item);
    return items;
}

QPixmap DrawFragmentSceneTool::createDefaultCursorPixmap() const
{
    HintPixmapFragmentItem frag_item(m_hint_item);
    return cursor_hint_from_graphics_item(&frag_item);
}

void DrawFragmentSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag actions are handled in onDragMove
        return;
    }
    auto [new_conf, overlay_atom] =
        getFragConfAndCoreAtomForScenePos(event->scenePos());
    m_hint_item.updateConformer(new_conf);
    std::unordered_set<RDKit::Bond*> bonds_to_mutate;
    if (overlay_atom != nullptr) {
        bonds_to_mutate = determine_fragment_bonds_to_mutate_to_single_bonds(
            m_frag, new_conf, *m_mol_model->getMol(), overlay_atom);
    }
    m_hint_item.updateSingleBondMutations(bonds_to_mutate);
    bool show_hint_structure = overlay_atom != nullptr;
    m_hint_item.setVisible(show_hint_structure);
    if (show_hint_structure == m_cursor_pixmap_shown) {
        // we want to hide the cursor pixmap while the full-sized hint structure
        // is shown
        emit newCursorHintRequested(
            show_hint_structure ? QPixmap() : getDefaultCursorPixmap());
        m_cursor_pixmap_shown = !show_hint_structure;
    }
}

void DrawFragmentSceneTool::onMouseLeave()
{
    StandardSceneToolBase::onMouseLeave();
    m_hint_item.hide();
}

void DrawFragmentSceneTool::onLeftButtonPress(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonPress(event);
    m_hint_item.show();
    // hide the cursor pixmap while we're showing the full-sized hint structure
    if (m_cursor_pixmap_shown) {
        emit newCursorHintRequested({});
        m_cursor_pixmap_shown = false;
    }
}

void DrawFragmentSceneTool::onLeftButtonRelease(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonRelease(event);
    emit newCursorHintRequested(getDefaultCursorPixmap());
    m_cursor_pixmap_shown = true;
}

void DrawFragmentSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonClick(event);
    auto [new_conf, overlay_atom] =
        getFragConfAndCoreAtomForScenePos(event->scenePos());
    addFragToModel(new_conf, overlay_atom);
}

void DrawFragmentSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragMove(event);
    auto maybe_conf = getConformerForDragToScenePos(event->scenePos());
    if (maybe_conf.has_value()) {
        // we only perform a drag action if the mouse press was over empty space
        m_hint_item.updateConformer(*maybe_conf);
    }
    m_hint_item.setVisible(maybe_conf.has_value());
}

void DrawFragmentSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragRelease(event);
    auto maybe_conf = getConformerForDragToScenePos(event->scenePos());
    if (maybe_conf.has_value()) {
        // we only perform a drag action if the mouse press was over empty space
        addFragToModel(*maybe_conf);
        m_hint_item.setVisible(false);
    }
}

std::pair<RDKit::Conformer, const RDKit::Atom*>
DrawFragmentSceneTool::getFragConfAndCoreAtomForScenePos(
    const QPointF& scene_pos) const
{
    auto* item = m_scene->getTopInteractiveItemAt(
        scene_pos, InteractiveItemFlag::MOLECULAR_NOT_AP);
    if (auto atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        const auto* core_atom = atom_item->getAtom();
        const auto& core = core_atom->getOwningMol();
        // if the atom has three or more neighbors, then don't try to add to it.
        // Instead, allow the fragment to follow the mouse cursor.
        if (core.getAtomDegree(core_atom) < 3) {
            auto new_conf = align_fragment_with_atom(m_frag, core_atom);
            return {new_conf, core_atom};
        }
    } else if (auto bond_item = qgraphicsitem_cast<BondItem*>(item)) {
        return align_fragment_with_bond(m_frag, bond_item->getBond());
    }
    // we're not over any part of the molecule, or we're over an atom with 3 or
    // more bonds
    auto mol_pos = to_mol_xy(scene_pos);
    auto new_conf = translate_fragment_center_to(m_frag, mol_pos);
    return {new_conf, nullptr};
}

std::optional<RDKit::Conformer>
DrawFragmentSceneTool::getConformerForDragToScenePos(
    const QPointF& scene_pos) const
{
    std::optional<RDKit::Conformer> maybe_conf;
    auto* item = m_scene->getTopInteractiveItemAt(
        m_mouse_press_scene_pos, InteractiveItemFlag::MOLECULAR_NOT_AP);
    if (item == nullptr) {
        // the drag was started from empty space, so we rotate the conformer
        auto mouse_press_mol_pos = to_mol_xy(m_mouse_press_scene_pos);
        // first, translate the conformer to where it was when the drag started
        maybe_conf = translate_fragment_center_to(m_frag, mouse_press_mol_pos);

        // find the centroid of the fragment, but ignore the attachment point,
        // which isn't being painted
        auto all_atoms = m_frag.atoms();
        std::unordered_set<const RDKit::Atom*> all_atoms_except_ap;
        std::copy_if(
            all_atoms.begin(), all_atoms.end(),
            std::inserter(all_atoms_except_ap, all_atoms_except_ap.end()),
            std::not_fn(is_attachment_point));
        auto centroid = find_centroid(*maybe_conf, all_atoms_except_ap);

        // determine the angle of the mouse cursor relative to the fragment
        // centroid, *not* where the drag started
        auto scene_centroid = to_scene_xy(centroid);
        auto new_angle = get_rounded_angle_radians(scene_centroid, scene_pos);
        rotate_conformer_radians(new_angle, centroid, *maybe_conf);
    }
    return maybe_conf;
}

void DrawFragmentSceneTool::addFragToModel(
    const RDKit::Conformer& conf, const RDKit::Atom* const overlay_atom)
{
    m_hint_item.hide();
    auto frag_copy = RDKit::ROMol(m_frag);
    frag_copy.getConformer() = conf;
    m_mol_model->addFragment(frag_copy, overlay_atom);
}

} // namespace sketcher
} // namespace schrodinger
