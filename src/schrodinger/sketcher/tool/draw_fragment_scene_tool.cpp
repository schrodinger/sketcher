#include "schrodinger/sketcher/tool/draw_fragment_scene_tool.h"

#include <stdexcept>

#include <boost/shared_ptr.hpp>

#include <QPointF>

#include <GraphMol/RWMol.h>

#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/fragment.h"
#include "schrodinger/sketcher/rdkit/molops.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

std::shared_ptr<DrawFragmentSceneTool>
get_draw_fragment_scene_tool(const RingTool& ring_tool, const Fonts& fonts,
                             const AtomItemSettings& atom_item_settings,
                             const BondItemSettings& bond_item_settings,
                             Scene* scene, MolModel* mol_model)
{
    std::string smiles = get_smiles_for_ring_tool(ring_tool);
    return get_draw_fragment_scene_tool(smiles, fonts, atom_item_settings,
                                        bond_item_settings, scene, mol_model);
}

std::shared_ptr<DrawFragmentSceneTool>
get_draw_fragment_scene_tool(const std::string& text_mol, const Fonts& fonts,
                             const AtomItemSettings& atom_item_settings,
                             const BondItemSettings& bond_item_settings,
                             Scene* scene, MolModel* mol_model)
{
    auto mol = text_to_mol(text_mol);
    return std::make_shared<DrawFragmentSceneTool>(
        *mol, fonts, atom_item_settings, bond_item_settings, scene, mol_model);
}

std::string get_smiles_for_ring_tool(const RingTool& ring_tool)
{
    switch (ring_tool) {
        case (RingTool::CYCLOPROPANE):
            return "*C1CC1 |$_AP1;;;$|";
        case (RingTool::CYCLOBUTANE):
            return "*C1CCC1 |$_AP1;;;;$|";
        case (RingTool::CYCLOPENTANE):
            return "*C1CCCC1 |$_AP1;;;;;$|";
        case (RingTool::CYCLOPENTADIENE):
            return "*C1=CCC=C1 |$_AP1;;;;;$|";
        case (RingTool::CYCLOHEXANE):
            return "*C1CCCCC1 |$_AP1;;;;;;$|";
        case (RingTool::BENZENE):
            return "*C1=CC=CC=C1 |$_AP1;;;;;;$|";
        case (RingTool::CYCLOHEPTANE):
            return "*C1CCCCCC1 |$_AP1;;;;;;;$|";
        case (RingTool::CYCLOOCTANE):
            return "*C1CCCCCCC1 |$_AP1;;;;;;;;$|";
    }
    throw std::runtime_error("Unrecognized ring type");
}

HintFragmentItem::HintFragmentItem(const RDKit::ROMol& fragment,
                                   const Fonts& fonts,
                                   const AtomItemSettings& atom_item_settings,
                                   const BondItemSettings& bond_item_settings,
                                   QGraphicsItem* parent) :
    QGraphicsItemGroup(parent),
    m_frag(fragment),
    m_atom_item_settings(atom_item_settings),
    m_bond_item_settings(bond_item_settings)
{
    m_atom_item_settings.setMonochromeColorScheme(HINT_COLOR);
    m_bond_item_settings.m_color = HINT_COLOR;
    m_bond_item_settings.m_bond_width = FRAGMENT_HINT_BOND_WIDTH;

    // We'll set this to visible once the mouse is over the scene
    setVisible(false);
    setZValue(static_cast<qreal>(ZOrder::HINT));

    auto [all_items, atom_to_atom_item, bond_to_bond_item] =
        create_graphics_items_for_mol(&m_frag, fonts, m_atom_item_settings,
                                      m_bond_item_settings,
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
}

void HintFragmentItem::updateConformer(const RDKit::Conformer& conformer)
{
    m_frag.getConformer() = conformer;
    update_conf_for_mol_graphics_items(m_atom_items, m_bond_items, m_frag);
}

DrawFragmentSceneTool::DrawFragmentSceneTool(
    const RDKit::ROMol& fragment, const Fonts& fonts,
    const AtomItemSettings& atom_item_settings,
    const BondItemSettings& bond_item_settings, Scene* scene,
    MolModel* mol_model) :
    AbstractSceneTool(scene, mol_model),
    m_frag(fragment),
    m_hint_item(m_frag, fonts, atom_item_settings, bond_item_settings)
{
}

std::vector<QGraphicsItem*> DrawFragmentSceneTool::getGraphicsItems()
{
    return {&m_hint_item};
}

void DrawFragmentSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    AbstractSceneTool::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag actions are handled in onDragMove
        return;
    }
    auto [new_conf, overlay_atom] =
        getFragConfAndCoreAtomForScenePos(event->scenePos());
    m_hint_item.updateConformer(new_conf);
    m_hint_item.show();
}

void DrawFragmentSceneTool::onMouseLeave()
{
    AbstractSceneTool::onMouseLeave();
    m_hint_item.hide();
}

void DrawFragmentSceneTool::onMouseClick(QGraphicsSceneMouseEvent* const event)
{
    AbstractSceneTool::onMouseClick(event);
    auto [new_conf, overlay_atom] =
        getFragConfAndCoreAtomForScenePos(event->scenePos());
    addFragToModel(new_conf, overlay_atom);
}

void DrawFragmentSceneTool::onDragMove(QGraphicsSceneMouseEvent* const event)
{
    AbstractSceneTool::onDragMove(event);
    auto maybe_conf = getConformerForDragToScenePos(event->scenePos());
    if (maybe_conf.has_value()) {
        // we only perform a drag action if the mouse press was over empty space
        m_hint_item.updateConformer(*maybe_conf);
    }
    m_hint_item.setVisible(maybe_conf.has_value());
}

void DrawFragmentSceneTool::onDragRelease(QGraphicsSceneMouseEvent* const event)
{
    AbstractSceneTool::onDragRelease(event);
    auto maybe_conf = getConformerForDragToScenePos(event->scenePos());
    if (maybe_conf.has_value()) {
        // we only perform a drag action if the mouse press was over empty space
        addFragToModel(*maybe_conf);
    }
}

std::pair<RDKit::Conformer, const RDKit::Atom*>
DrawFragmentSceneTool::getFragConfAndCoreAtomForScenePos(
    const QPointF& scene_pos) const
{
    auto* item = m_scene->getTopInteractiveItemAt(
        scene_pos, InteractiveItemFlag::MOLECULAR_NOT_AP);
    if (auto atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
        const auto core_atom = atom_item->getAtom();
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
    auto new_conf = translate_fragment(m_frag, mol_pos);
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
        maybe_conf = translate_fragment(m_frag, mouse_press_mol_pos);

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
