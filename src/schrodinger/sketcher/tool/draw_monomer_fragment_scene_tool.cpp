#include "schrodinger/sketcher/tool/draw_monomer_fragment_scene_tool.h"

#include <ranges>

#include <fmt/core.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/monomer_directions.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/monomer_hint_fragment_item.h"
#include "schrodinger/sketcher/molviewer/monomer_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"

namespace schrodinger
{
namespace sketcher
{

// The HELM string we use to generate the nucleotide fragments and the
// associated atom indices of the monomers.
constexpr std::string_view NT_HELM_FMT = "RNA1{{{}.{}({})}}$$$$V2.0";
constexpr unsigned int PHOS_IDX = 0;
constexpr unsigned int SUGAR_IDX = 1;
constexpr unsigned int BASE_IDX = 2;

/**
 * We round the drag angle to 8 increments, i.e., 45 degree increments, instead
 * of the 6 increments used for atomistic drags.
 */
constexpr int DRAG_ANGLE_ROUNDING = 8;

/**
 * The allowable interactions used for the nucleotide fragments. We can add a
 * nucleotide to the 5' side of a phosphate, the 3' side of a sugar, or paired
 * with a base.
 */
const std::vector<MonomerFragmentAttachmentInfo> NUCLEOTIDE_ATTACHMENT_INFO = {
    {MonomerType::NA_PHOSPHATE, ap_model_name_for(NAPhosphateAP::TO_PREV_SUGAR),
     SUGAR_IDX, ap_model_name_for(NASugarAP::THREE_PRIME)},
    {MonomerType::NA_SUGAR, ap_model_name_for(NASugarAP::THREE_PRIME), PHOS_IDX,
     ap_model_name_for(NAPhosphateAP::TO_PREV_SUGAR)},
    {MonomerType::NA_BASE, NA_BASE_AP_PAIR, BASE_IDX, NA_BASE_AP_PAIR}};

/**
 * The rotation angle adjustment required so that the base monomer follows the
 * cursor during a click-and-drag of the nucleotide scene tool
 */
constexpr double NUCLEOTIDE_DRAG_ANGLE_ADJUSTMENT = M_PI / 2.0;

std::shared_ptr<DrawMonomerFragmentSceneTool>
get_nucleotide_fragment_scene_tool(const std::string& sugar,
                                   const std::string& base,
                                   const std::string& phos, const Fonts& fonts,
                                   Scene* scene, MolModel* mol_model)
{
    auto to_helm_monomer = [](const std::string& monomer) {
        if (monomer.size() > 1) {
            return "[" + monomer + "]";
        } else {
            return monomer;
        }
    };

    auto nt_helm = fmt::format(NT_HELM_FMT, to_helm_monomer(phos),
                               to_helm_monomer(sugar), to_helm_monomer(base));
    auto mol =
        rdkit_extensions::to_rdkit(nt_helm, rdkit_extensions::Format::HELM);
    prepare_mol(*mol);
    return std::make_shared<DrawMonomerFragmentSceneTool>(
        *mol, NUCLEOTIDE_ATTACHMENT_INFO, SUGAR_IDX,
        NUCLEOTIDE_DRAG_ANGLE_ADJUSTMENT, fonts, scene, mol_model);
}

/**
 * Modify the given molecule in place to translate the atom with index
 * atom_idx_to_center_on_coords to coords, and then rotate it by the specified
 * angle.
 */
static void move_mol_to_coords_and_rotate(
    RDKit::ROMol& mol, const unsigned int atom_idx_to_center_on_coords,
    const RDGeom::Point3D coords, const double angle_radians)
{
    auto& conf = mol.getConformer();
    auto mol_pos = conf.getAtomPos(atom_idx_to_center_on_coords);
    auto coord_adjustment = coords - mol_pos;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        conf.setAtomPos(i, conf.getAtomPos(i) + coord_adjustment);
    }
    if (angle_radians != 0.0) {
        rotate_conformer_radians(angle_radians, coords, conf);
    }
}

DrawMonomerFragmentSceneTool::DrawMonomerFragmentSceneTool(
    const RDKit::ROMol& mol,
    const std::vector<MonomerFragmentAttachmentInfo> attachment_info,
    const int index_to_center_on_click, const double drag_angle_adjustment,
    const Fonts& fonts, Scene* scene, MolModel* mol_model) :
    AbstractMonomerSceneTool(fonts, scene, mol_model),
    m_frag(mol),
    m_attachment_info(attachment_info),
    m_index_to_center_on_click(index_to_center_on_click),
    m_drag_angle_adjustment(drag_angle_adjustment)
{
    set_all_atoms_monomeric(m_frag);
    m_frag.getConformer().set3D(false);
}

void DrawMonomerFragmentSceneTool::onMouseMove(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag logic is handled in onDragMove
        return;
    }
    QPointF scene_pos = event->scenePos();
    auto hovered_monomer_item = getTopMonomerItemAt(scene_pos);
    QGraphicsItem* hovered_item = hovered_monomer_item;
    if (hovered_monomer_item == nullptr) {
        // if we're not over a monomer, check to see if we're over a connector
        hovered_item = getTopMonomerConnectorItemAt(scene_pos);
    }

    if (hovered_item != m_hovered_item) {
        // we're hovering over something new, so update the attachment point
        // labels
        m_hovered_item = hovered_item;
        drawAttachmentPointLabelsFor(hovered_item);
    }

    auto* hovered_ap_item =
        getConnectableUnboundAttachmentPointAt(scene_pos, hovered_monomer_item);
    if (hovered_ap_item != m_hovered_ap_item) {
        // update which attachment point is hovered
        m_hovered_ap_item = hovered_ap_item;
        // hide the hovered "nubbin" and draw a hint fragment instead
        for (auto* ap_item : m_unbound_ap_items) {
            ap_item->setVisible(hovered_ap_item == nullptr ||
                                ap_item != hovered_ap_item);
        }
        drawFragmentHintFor(hovered_monomer_item, hovered_ap_item);
        if (hovered_ap_item != nullptr) {
            // only highlight the monomer if a click would actually build the
            // fragment
            m_predictive_highlighting_item.highlightItem(hovered_monomer_item);
        }
    }

    if (hovered_ap_item == nullptr) {
        m_predictive_highlighting_item.clearHighlightingPath();
    }

    // we want to show either the hint fragment or the cursor hint, but not both
    bool drew_hint_fragment = hovered_ap_item != nullptr;
    setCursorHintShown(!drew_hint_fragment);
}

void DrawMonomerFragmentSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonClick(event);
    QPointF scene_pos = event->scenePos();
    auto hovered_monomer_item = getTopMonomerItemAt(scene_pos);
    auto* hovered_ap_item =
        getConnectableUnboundAttachmentPointAt(scene_pos, hovered_monomer_item);
    if (hovered_ap_item != nullptr) {
        // the user clicked on a valid attachment point (or a monomer with a
        // valid attachment point), so add the fragment and connect it to the
        // specified attachment point
        addBoundFragmentToMolModel(hovered_monomer_item, hovered_ap_item);
    } else if (hovered_monomer_item == nullptr) {
        // the user clicked on empty space, so add the fragment in its default
        // configuration
        addUnboundFragmentToMolModel(scene_pos);
    }
}

void DrawMonomerFragmentSceneTool::onLeftButtonPress(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonPress(event);
    auto scene_pos = event->scenePos();
    auto hovered_monomer = getTopMonomerItemAt(scene_pos);
    m_press_was_over_monomer = hovered_monomer != nullptr;
    if (!m_press_was_over_monomer) {
        setCursorHintShown(false);

        // create the structure hint
        auto frag_copy = std::make_shared<RDKit::RWMol>(m_frag);
        move_mol_to_coords_and_rotate(*frag_copy, m_index_to_center_on_click,
                                      to_mol_xy(scene_pos), 0.0);
        m_hint_fragment_item = new MonomerHintFragmentItem(
            frag_copy, m_fonts, {}, -1, m_monomer_background_color);
        m_scene->addItem(m_hint_fragment_item);
    }
}

void DrawMonomerFragmentSceneTool::onLeftButtonRelease(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonRelease(event);
    setCursorHintShown(true);
    delete m_hint_fragment_item;
}

void DrawMonomerFragmentSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonDragMove(event);
    if (m_press_was_over_monomer) {
        return;
    }
    auto drag_angle = getDragAngle(event->scenePos());
    m_hint_fragment_item->setRotation(drag_angle, m_index_to_center_on_click);
}

void DrawMonomerFragmentSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonDragRelease(event);
    if (m_press_was_over_monomer) {
        return;
    }
    auto drag_angle = getDragAngle(event->scenePos());
    addUnboundFragmentToMolModel(m_mouse_press_scene_pos, drag_angle);
}

UnboundMonomericAttachmentPointItem*
DrawMonomerFragmentSceneTool::getConnectableUnboundAttachmentPointAt(
    const QPointF& scene_pos, const AbstractMonomerItem* hovered_monomer) const
{
    if (m_unbound_ap_items.empty() || hovered_monomer == nullptr) {
        return nullptr;
    }

    auto hovered_monomer_type = get_monomer_type(hovered_monomer->getAtom());
    // for this monomer type, figure out which attachment point we can handle
    // (i.e which attachment points appear in m_attachment_info)
    std::vector<std::string> acceptable_ap_model_names;
    for (auto attachment_info : m_attachment_info) {
        if (attachment_info.struc_monomer_type == hovered_monomer_type) {
            acceptable_ap_model_names.push_back(
                attachment_info.struc_ap_model_name);
        }
    }

    // if the user is hovering over an attachment point, figure out whether it's
    // one we can handle
    for (auto* ap_item : m_unbound_ap_items) {
        auto ap_model_name = ap_item->getAttachmentPoint().model_name;
        if (ap_item->withinHoverArea(scene_pos)) {
            if (std::ranges::count(acceptable_ap_model_names, ap_model_name)) {
                // the user is hovering over an attachment point that appears in
                // our attachment info, so we can connect to it
                return ap_item;
            } else {
                // The user is hovering over an attachment point, but it doesn't
                // appear in our attachment info. We ignore it since we can't
                // connect to it.
                return nullptr;
            }
        }
    }

    // the user is hovering over the monomer itself, not an attachment point, so
    // find the attachment point with the lowest index in
    // acceptable_ap_model_names (i.e. out of all the attachment points
    // available on this monomer, which one (if any) appears first in
    // m_attachment_info)
    ptrdiff_t best_dist = acceptable_ap_model_names.size();
    UnboundMonomericAttachmentPointItem* best_ap_item = nullptr;
    for (auto* ap_item : m_unbound_ap_items) {
        auto it = std::ranges::find(acceptable_ap_model_names,
                                    ap_item->getAttachmentPoint().model_name);
        auto it_dist = std::distance(acceptable_ap_model_names.begin(), it);
        if (it_dist < best_dist) {
            best_dist = it_dist;
            best_ap_item = ap_item;
        }
    }
    // this will be a nullptr if none of the available attachment points were
    // acceptable
    return best_ap_item;
}

/**
 * @return the attachment info entry corresponding to the type of the given
 * monomer and the specified attachment point
 * @throw runtime_error if no such attachment info is found
 */
static MonomerFragmentAttachmentInfo
get_attachment_info(std::vector<MonomerFragmentAttachmentInfo> attachment_info,
                    const RDKit::Atom* const monomer,
                    const UnboundAttachmentPoint& ap)
{
    auto monomer_type = get_monomer_type(monomer);

    auto attachment_info_it = std::ranges::find_if(
        attachment_info, [monomer_type, &ap](auto attachment_info) {
            return attachment_info.struc_monomer_type == monomer_type &&
                   attachment_info.struc_ap_model_name == ap.model_name;
        });
    if (attachment_info_it == attachment_info.end()) {
        throw std::runtime_error(
            "Attachment point not valid for this fragment");
    }
    return *attachment_info_it;
}

void DrawMonomerFragmentSceneTool::drawFragmentHintFor(
    const AbstractMonomerItem* hovered_item,
    const UnboundMonomericAttachmentPointItem* hovered_ap_item)
{
    clearHintFragmentItem();
    if (hovered_ap_item == nullptr) {
        return;
    }

    auto* hovered_monomer = hovered_item->getAtom();
    auto hovered_ap = hovered_ap_item->getAttachmentPoint();
    auto attachment_info =
        get_attachment_info(m_attachment_info, hovered_monomer, hovered_ap);
    auto scene_struc_conf = hovered_monomer->getOwningMol().getConformer();
    auto hovered_monomer_pos =
        scene_struc_conf.getAtomPos(hovered_monomer->getIdx());
    auto hint_mol =
        getPositionedFragment(hovered_monomer, hovered_ap, attachment_info);

    // Add a copy of the hovered monomer to the hint fragment along with the
    // appropriate connection. (We'll hide the graphics item for this monomer,
    // but we need it so we can create the connection.)
    auto hovered_monomer_copy_idx =
        hint_mol->addAtom(new RDKit::Atom(*hovered_monomer), false, true);
    hint_mol->getConformer().setAtomPos(hovered_monomer_copy_idx,
                                        hovered_monomer_pos);
    auto new_connection_idx = add_monomer_connection(
        *hint_mol, hovered_monomer_copy_idx,
        attachment_info.struc_ap_model_name, attachment_info.frag_monomer_idx,
        attachment_info.frag_monomer_ap_model_name);

    m_hint_fragment_item = new MonomerHintFragmentItem(
        hint_mol, m_fonts, {hovered_monomer_copy_idx}, new_connection_idx,
        m_monomer_background_color);
    m_scene->addItem(m_hint_fragment_item);
}

/**
 * @return the attachment point with the specified model name
 * @throw runtime_error if no such attachment point is found
 */
static UnboundAttachmentPoint find_unbound_attachment_point(
    const std::vector<UnboundAttachmentPoint>& unbound_attachment_points,
    const std::string_view ap_model_name)
{
    auto unbound_ap_it = std::ranges::find_if(
        unbound_attachment_points,
        [&ap_model_name](auto ap) { return ap.model_name == ap_model_name; });
    if (unbound_ap_it == unbound_attachment_points.end()) {
        throw std::runtime_error(
            "Attachment point not valid for this fragment");
    }
    return *unbound_ap_it;
}

std::shared_ptr<RDKit::RWMol>
DrawMonomerFragmentSceneTool::getPositionedFragment(
    const RDKit::Atom* const hovered_monomer,
    const UnboundAttachmentPoint& hovered_ap,
    const MonomerFragmentAttachmentInfo& attachment_info)
{
    // figure out where to translate the fragment to
    const auto& conf = hovered_monomer->getOwningMol().getConformer();
    auto monomer_pos = conf.getAtomPos(hovered_monomer->getIdx());
    auto hint_pos =
        get_default_coords_for_bound_monomer(monomer_pos, hovered_ap.direction);

    // figure out what angle we need to rotate the fragment so that the
    // attachment points to be connected will line up
    const auto* bound_frag_monomer =
        m_frag.getAtomWithIdx(attachment_info.frag_monomer_idx);
    auto [bound_attachment_points, unbound_attachment_points] =
        get_attachment_points_for_monomer(bound_frag_monomer);
    auto unbound_ap = find_unbound_attachment_point(
        unbound_attachment_points, attachment_info.frag_monomer_ap_model_name);
    auto vec_from_frag_ap =
        rdkit_extensions::direction_to_vector(unbound_ap.direction);
    auto vec_to_hovered_ap =
        -rdkit_extensions::direction_to_vector(hovered_ap.direction);
    auto angle_to_rotate =
        get_angle_radians(vec_to_hovered_ap, {0, 0, 0}, vec_from_frag_ap);

    auto frag_copy = std::make_shared<RDKit::RWMol>(m_frag);
    move_mol_to_coords_and_rotate(*frag_copy, attachment_info.frag_monomer_idx,
                                  hint_pos, angle_to_rotate);
    return frag_copy;
}

void DrawMonomerFragmentSceneTool::addBoundFragmentToMolModel(
    const AbstractMonomerItem* hovered_monomer_item,
    const UnboundMonomericAttachmentPointItem* hovered_ap_item)
{
    auto* hovered_monomer = hovered_monomer_item->getAtom();
    auto hovered_monomer_idx = hovered_monomer->getIdx();
    auto hovered_ap = hovered_ap_item->getAttachmentPoint();
    auto attachment_info =
        get_attachment_info(m_attachment_info, hovered_monomer, hovered_ap);
    auto positioned_frag =
        getPositionedFragment(hovered_monomer, hovered_ap, attachment_info);

    auto first_frag_idx = m_mol_model->getMol()->getNumAtoms();
    auto undo_raii = m_mol_model->createUndoMacro("Add monomeric fragment");
    m_mol_model->addMol(*positioned_frag, "Add monomeric fragment",
                        /* reposition_mol = */ false,
                        /* new_mol_added = */ false,
                        /* enforce_size_limit = */ false);
    auto* mol = m_mol_model->getMol();
    hovered_monomer = mol->getAtomWithIdx(hovered_monomer_idx);
    auto frag_monomer_to_bind_idx =
        first_frag_idx + attachment_info.frag_monomer_idx;
    auto frag_monomer_to_bind = mol->getAtomWithIdx(frag_monomer_to_bind_idx);
    m_mol_model->addMonomericConnection(
        hovered_monomer, hovered_ap.model_name, frag_monomer_to_bind,
        attachment_info.frag_monomer_ap_model_name);
}

void DrawMonomerFragmentSceneTool::addUnboundFragmentToMolModel(
    const QPointF& scene_pos, const double rotation)
{
    auto frag_copy = std::make_shared<RDKit::RWMol>(m_frag);
    move_mol_to_coords_and_rotate(*frag_copy, m_index_to_center_on_click,
                                  to_mol_xy(scene_pos), rotation);
    m_mol_model->addMol(*frag_copy, "Add monomeric fragment",
                        /* reposition_mol = */ false,
                        /* new_mol_added = */ false,
                        /* enforce_size_limit = */ false);
}

double DrawMonomerFragmentSceneTool::getDragAngle(const QPointF& scene_pos)
{
    return get_rounded_angle_radians(m_mouse_press_scene_pos, scene_pos,
                                     DRAG_ANGLE_ROUNDING) +
           m_drag_angle_adjustment;
}

} // namespace sketcher
} // namespace schrodinger
