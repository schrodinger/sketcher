#include "schrodinger/sketcher/tool/abstract_draw_monomer_or_monomeric_connection_scene_tool.h"

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

#include <QtAssert>
#include <QtMath>
#include <QGraphicsItem>
#include <QGraphicsSimpleTextItem>
#include <QPainterPath>
#include <QPen>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/monomer_directions.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/monomer_attachment_point_labels.h"
#include "schrodinger/sketcher/molviewer/monomer_connector_item.h"
#include "schrodinger/sketcher/molviewer/monomer_hint_fragment_item.h"
#include "schrodinger/sketcher/molviewer/monomer_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/unbound_monomeric_attachment_point_item.h"

namespace schrodinger
{
namespace sketcher
{

using rdkit_extensions::Direction;

AbstractDrawMonomerOrMonomericConnectionSceneTool::
    AbstractDrawMonomerOrMonomericConnectionSceneTool(
        const std::string& res_name,
        const rdkit_extensions::ChainType chain_type, const Fonts& fonts,
        Scene* scene, MolModel* mol_model) :
    AbstractMonomerSceneTool(fonts, scene, mol_model),
    m_res_name(res_name),
    m_chain_type(chain_type),
    m_bolded_fonts(fonts)
{
    // we handle predictive highlighting manually in onMouseMove, so we disable
    // StandardSceneToolsBase's predictive highlighting
    m_highlight_types = InteractiveItemFlag::NONE;
    // make sure that the cursor hint font is more easily readable at small
    // size. (Note that m_bolded_fonts is a copy of the Scene's fonts, so this
    // change won't affect anything else.)
    m_bolded_fonts.m_main_label_font.setBold(true);
    m_bolded_fonts.updateFontMetrics();
    if (chain_type == rdkit_extensions::ChainType::PEPTIDE) {
        m_monomer_type = MonomerType::PEPTIDE;
    } else if (chain_type == rdkit_extensions::ChainType::CHEM) {
        m_monomer_type = MonomerType::CHEM;
    } else {
        m_monomer_type = get_na_monomer_type_from_res_name(res_name);
    }
}

AbstractDrawMonomerOrMonomericConnectionSceneTool::
    ~AbstractDrawMonomerOrMonomericConnectionSceneTool()
{
    // explicitly erase any drag end attachment point labels. Without this, the
    // bound attachment point labels would be deleted regardless when
    // m_drag_end_attachment_point_labels_group is destroyed, but the unbound
    // attachment point labels are parented to their monomer, so they would
    // outlive the scene tool without this call.
    clearDragEndAttachmentPointsLabels();
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::onStructureUpdated()
{
    AbstractMonomerSceneTool::onStructureUpdated();
    clear_graphics_item_group(m_drag_end_attachment_point_labels_group);
    // any drag-end attachment point graphics items have already been deleted in
    // response to their parent monomer graphics item being deleted due to the
    // structure update. We're just clear the stale pointers here
    m_drag_end_unbound_ap_items.clear();
    m_drag_start_monomer_item = nullptr;
    m_drag_end_monomer_item = nullptr;
}

std::vector<QGraphicsItem*>
AbstractDrawMonomerOrMonomericConnectionSceneTool::getGraphicsItems()
{
    auto items = AbstractMonomerSceneTool::getGraphicsItems();
    items.push_back(&m_drag_end_attachment_point_labels_group);
    return items;
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::
    updateColorsAfterBackgroundColorChange(bool is_dark_mode)
{
    AbstractMonomerSceneTool::updateColorsAfterBackgroundColorChange(
        is_dark_mode);
    m_monomer_background_color =
        is_dark_mode ? DARK_BACKGROUND_COLOR : LIGHT_BACKGROUND_COLOR;
    m_drag_end_inactive_ap_color =
        is_dark_mode ? DRAG_END_INACTIVE_AP_DARK_BG : DRAG_END_INACTIVE_AP;
}

UnboundMonomericAttachmentPointItem* find_preferred_attachment_point_by_num(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items,
    const std::vector<int>& preferred_order)
{
    auto index_of = [&preferred_order](const auto* ap_item) {
        auto it = std::ranges::find(preferred_order,
                                    ap_item->getAttachmentPoint().num);
        auto dist = std::distance(preferred_order.begin(), it);
        // we know that dist is positive
        return static_cast<std::size_t>(dist);
    };
    auto min_it = std::ranges::min_element(
        unbound_ap_items,
        [&index_of](const auto* ap_item_left, const auto* ap_item_right) {
            return index_of(ap_item_left) < index_of(ap_item_right);
        });
    // make sure that the attachment point we found is actually on the
    // preferred_order list
    if (index_of(*min_it) < preferred_order.size()) {
        return *min_it;
    }
    return nullptr;
}

UnboundMonomericAttachmentPointItem* find_min_attachment_point_by_num(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items)
{
    auto min_it = std::ranges::min_element(
        unbound_ap_items,
        [](const auto* ap_item_left, const auto* ap_item_right) {
            auto left_num = ap_item_left->getAttachmentPoint().num;
            auto right_num = ap_item_right->getAttachmentPoint().num;
            return right_num < 0 || (left_num >= 0 && left_num < right_num);
        });
    return *min_it;
}

UnboundMonomericAttachmentPointItem* find_attachment_point_with_name(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items,
    const std::string_view name)
{
    auto it =
        std::ranges::find_if(unbound_ap_items, [&name](const auto* ap_item) {
            return (ap_item->getAttachmentPoint().model_name == name);
        });
    if (it == unbound_ap_items.end()) {
        return nullptr;
    }
    return *it;
}

/**
 * @return the unbound attachment point graphics item representing an attachment
 * point with the given number. Will return nullptr if no such attachment point
 * is found.
 */
[[nodiscard]] static UnboundMonomericAttachmentPointItem*
find_attachment_point_with_num(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items,
    const int num)
{
    auto it =
        std::ranges::find_if(unbound_ap_items, [&num](const auto* ap_item) {
            return (ap_item->getAttachmentPoint().num == num);
        });
    if (it == unbound_ap_items.end()) {
        return nullptr;
    }
    return *it;
}

UnboundMonomericAttachmentPointItem* get_default_attachment_point(
    const MonomerType hovered_type, const MonomerType tool_type,
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items)
{
    if (unbound_ap_items.empty()) {
        return nullptr;
    } else if (hovered_type == MonomerType::CHEM) {
        return find_min_attachment_point_by_num(unbound_ap_items);
    } else if (hovered_type == MonomerType::PEPTIDE) {
        if (tool_type == MonomerType::PEPTIDE) {
            return find_preferred_attachment_point_by_num(
                unbound_ap_items,
                {PeptideAP::C, PeptideAP::N, PeptideAP::SIDECHAIN});
        } else if (tool_type == MonomerType::CHEM) {
            return find_attachment_point_with_num(unbound_ap_items,
                                                  PeptideAP::SIDECHAIN);
        }
    } else if (hovered_type == MonomerType::NA_BASE) {
        if (tool_type == MonomerType::NA_BASE ||
            tool_type == MonomerType::CHEM) {
            return find_attachment_point_with_name(unbound_ap_items,
                                                   H_BOND_AP_MODEL_NAME);
        } else if (tool_type == MonomerType::NA_SUGAR) {
            return find_attachment_point_with_num(unbound_ap_items,
                                                  NA_BASE_AP_N1_9);
        }
    } else if (hovered_type == MonomerType::NA_SUGAR) {
        if (tool_type == MonomerType::NA_BASE) {
            return find_attachment_point_with_num(unbound_ap_items,
                                                  NASugarAP::ONE_PRIME);
        } else if (tool_type == MonomerType::NA_PHOSPHATE) {
            return find_preferred_attachment_point_by_num(
                unbound_ap_items,
                {NASugarAP::THREE_PRIME, NASugarAP::FIVE_PRIME});
        }
    } else if (hovered_type == MonomerType::NA_PHOSPHATE) {
        if (tool_type == MonomerType::NA_SUGAR) {
            return find_preferred_attachment_point_by_num(
                unbound_ap_items,
                {NAPhosphateAP::TO_NEXT_SUGAR, NAPhosphateAP::TO_PREV_SUGAR});
        }
    }
    return nullptr;
}

std::string
get_attachment_point_for_new_monomer(const MonomerType existing_monomer_type,
                                     const std::string_view existing_monomer_ap,
                                     const MonomerType new_monomer_type)
{
    switch (new_monomer_type) {
        case MonomerType::CHEM:
            return "R1";
        case MonomerType::PEPTIDE:
            if (existing_monomer_type == MonomerType::PEPTIDE) {
                if (existing_monomer_ap == ap_model_name_for(PeptideAP::N)) {
                    return ap_model_name_for(PeptideAP::C);
                } else if (existing_monomer_ap ==
                           ap_model_name_for(PeptideAP::C)) {
                    return ap_model_name_for(PeptideAP::N);
                }
            }
            return ap_model_name_for(PeptideAP::SIDECHAIN);
        case MonomerType::NA_BASE:
            if (existing_monomer_type == MonomerType::NA_SUGAR &&
                existing_monomer_ap ==
                    ap_model_name_for(NASugarAP::ONE_PRIME)) {
                return ap_model_name_for(NA_BASE_AP_N1_9);
            }
            return H_BOND_AP_MODEL_NAME;
        case MonomerType::NA_SUGAR:
            if (existing_monomer_type == MonomerType::NA_PHOSPHATE) {
                if (existing_monomer_ap ==
                    ap_model_name_for(NAPhosphateAP::TO_PREV_SUGAR)) {
                    return ap_model_name_for(NASugarAP::THREE_PRIME);
                } else if (existing_monomer_ap ==
                           ap_model_name_for(NAPhosphateAP::TO_NEXT_SUGAR)) {
                    return ap_model_name_for(NASugarAP::FIVE_PRIME);
                }
            }
            return ap_model_name_for(NASugarAP::ONE_PRIME);
        case MonomerType::NA_PHOSPHATE:
            if (existing_monomer_type == MonomerType::NA_SUGAR &&
                existing_monomer_ap ==
                    ap_model_name_for(NASugarAP::THREE_PRIME)) {
                return ap_model_name_for(NAPhosphateAP::TO_PREV_SUGAR);
            }
            return ap_model_name_for(NAPhosphateAP::TO_NEXT_SUGAR);
        default:
            Q_UNREACHABLE_RETURN("");
    }
}

UnboundMonomericAttachmentPointItem*
AbstractDrawMonomerOrMonomericConnectionSceneTool::getUnboundAttachmentPointAt(
    const QPointF& scene_pos,
    const bool no_default_if_click_should_mutate) const
{
    if (m_unbound_ap_items.empty()) {
        return nullptr;
    }
    for (auto* ap_item : m_unbound_ap_items) {
        if (ap_item->withinHoverArea(scene_pos)) {
            return ap_item;
        }
    }
    return getDefaultUnboundAttachmentPointForHoveredMonomer(
        no_default_if_click_should_mutate);
}

std::tuple<const RDKit::Atom*, MonomerType>
get_monomer_and_type(const AbstractMonomerItem* const monomer_item)
{
    auto* monomer = monomer_item->getAtom();
    auto monomer_type = get_monomer_type(monomer);
    return {monomer, monomer_type};
}

UnboundMonomericAttachmentPointItem*
AbstractDrawMonomerOrMonomericConnectionSceneTool::
    getUnboundDragEndAttachmentPointAt(const QPointF& scene_pos) const
{
    if (m_drag_end_unbound_ap_items.empty()) {
        return nullptr;
    }
    for (auto* ap_item : m_drag_end_unbound_ap_items) {
        if (ap_item->withinHoverArea(scene_pos)) {
            return ap_item;
        }
    }

    // get the monomer type for the start and end monomers
    auto get_drag_monomer_type =
        [this](const AbstractMonomerItem* const monomer_item) {
            if (monomer_item == nullptr) {
                return m_monomer_type;
            }
            auto [monomer, monomer_type] = get_monomer_and_type(monomer_item);
            return monomer_type;
        };
    auto drag_start_monomer_type =
        get_drag_monomer_type(m_drag_start_monomer_item);
    auto drag_end_monomer_type = get_drag_monomer_type(m_drag_end_monomer_item);

    // if the user is hovered over the monomer itself (not an attachment point)
    // and the "correct" attachment point is available (e.g. N terminus when we
    // dragged from a C terminus), use that one
    auto ideal_ap_model_name = get_attachment_point_for_new_monomer(
        drag_start_monomer_type, m_drag_start_ap_model_name,
        drag_end_monomer_type);
    for (auto* ap_item : m_drag_end_unbound_ap_items) {
        if (ap_item->getAttachmentPoint().model_name == ideal_ap_model_name) {
            return ap_item;
        }
    }

    // if the correct attachment point isn't available then there's probably no
    // good option, so use whatever we would've defaulted if we'd started the
    // drag at this monomer.
    return get_default_attachment_point(drag_end_monomer_type,
                                        drag_start_monomer_type,
                                        m_drag_end_unbound_ap_items);
}

std::tuple<const RDKit::Atom*, MonomerType>
AbstractDrawMonomerOrMonomericConnectionSceneTool::getHoveredMonomerAndType()
    const
{
    if (!item_matches_type_flag(m_hovered_monomeric_item,
                                InteractiveItemFlag::MONOMER)) {
        throw std::runtime_error("No hovered monomer");
    }
    const auto* monomer_item =
        static_cast<const AbstractMonomerItem*>(m_hovered_monomeric_item);
    return get_monomer_and_type(monomer_item);
}

UnboundMonomericAttachmentPointItem*
AbstractDrawMonomerOrMonomericConnectionSceneTool::
    getDefaultUnboundAttachmentPointForHoveredMonomer(
        const bool no_default_if_click_should_mutate) const
{
    auto [monomer, monomer_type] = getHoveredMonomerAndType();
    if (no_default_if_click_should_mutate &&
        clickShouldMutate(monomer, monomer_type)) {
        return nullptr;
    }
    return get_default_attachment_point(monomer_type, m_monomer_type,
                                        m_unbound_ap_items);
}

bool AbstractDrawMonomerOrMonomericConnectionSceneTool::
    shouldShowPredictiveHighlighting() const
{

    if (m_hovered_monomeric_item == nullptr ||
        !item_matches_type_flag(m_hovered_monomeric_item,
                                InteractiveItemFlag::MONOMER)) {
        return false;
    }
    auto [monomer, monomer_type] = getHoveredMonomerAndType();
    if (clickShouldMutate(monomer, monomer_type)) {
        return true;
    }
    return get_default_attachment_point(monomer_type, m_monomer_type,
                                        m_unbound_ap_items) != nullptr;
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::onMouseMove(
    QGraphicsSceneMouseEvent* const event)
{
    // the base class call updates m_hovered_monomeric_item and the attachment
    // point labels
    auto prev_hovered_monomeric_item = m_hovered_monomeric_item;
    AbstractMonomerSceneTool::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag logic is handled in onDragMove
        return;
    }
    QPointF scene_pos = event->scenePos();

    if (m_hovered_monomeric_item != prev_hovered_monomeric_item) {
        if (shouldShowPredictiveHighlighting()) {
            m_predictive_highlighting_item.highlightItem(
                m_hovered_monomeric_item);
        } else {
            m_predictive_highlighting_item.clearHighlightingPath();
        }
    }

    auto* hovered_ap_item = getUnboundAttachmentPointAt(scene_pos, true);
    if (hovered_ap_item != m_hovered_ap_item) {
        // update which attachment point is hovered
        m_hovered_ap_item = hovered_ap_item;
        updateHoveredUnboundAP(hovered_ap_item);
    }

    // we want to show either the hint fragment or the cursor hint, but not both
    bool drew_hint_fragment = hovered_ap_item != nullptr;
    setCursorHintShown(!drew_hint_fragment);
}

static RDGeom::Point3D get_coords_for_monomer(const RDKit::Atom* const monomer)
{
    auto& conf = monomer->getOwningMol().getConformer();
    return conf.getAtomPos(monomer->getIdx());
}

RDGeom::Point3D
get_default_coords_for_bound_monomer(const RDKit::Atom* const monomer,
                                     const Direction dir)
{
    auto monomer_pos = get_coords_for_monomer(monomer);
    return get_default_coords_for_bound_monomer(monomer_pos, dir);
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::createHintFragmentItem(
    HintFragmentMonomerInfo& monomer_one_info,
    HintFragmentMonomerInfo& monomer_two_info)
{
    auto frag = std::make_shared<RDKit::RWMol>();

    // create the two monomers
    auto monomer_one_idx =
        frag->addAtom(monomer_one_info.monomer.release(), true, true);
    auto monomer_two_idx =
        frag->addAtom(monomer_two_info.monomer.release(), true, true);
    set_all_atoms_monomeric(*frag);
    // create the connection between them
    auto bond_index_to_label = add_monomer_connection(
        *frag, monomer_one_idx, monomer_one_info.ap_model_name, monomer_two_idx,
        monomer_two_info.ap_model_name);

    // Add a conformer with the atom coordinates
    auto* frag_conf = new RDKit::Conformer(frag->getNumAtoms());
    frag_conf->set3D(false);
    frag_conf->setAtomPos(monomer_one_idx, monomer_one_info.pos);
    frag_conf->setAtomPos(monomer_two_idx, monomer_two_info.pos);
    frag->addConformer(frag_conf, true);

    // hide the monomers that already exist in the Scene
    std::vector<size_t> atom_indices_to_hide;
    if (monomer_one_info.atom_idx >= 0) {
        atom_indices_to_hide.push_back(monomer_one_idx);
    }
    if (monomer_two_info.atom_idx >= 0) {
        atom_indices_to_hide.push_back(monomer_two_idx);
    }

    m_hint_fragment_item = new MonomerHintFragmentItem(
        frag, m_bolded_fonts, atom_indices_to_hide, bond_index_to_label,
        m_monomer_background_color);
    m_scene->addItem(m_hint_fragment_item);
}

std::string AbstractDrawMonomerOrMonomericConnectionSceneTool::
    getDefaultDragStartAPModelName() const
{
    return m_monomer_type == MonomerType::NA_BASE
               ? H_BOND_AP_MODEL_NAME
               : ap_model_name_for(PeptideAP::C);
}

HintFragmentMonomerInfo AbstractDrawMonomerOrMonomericConnectionSceneTool::
    createHintFragmentMonomerInfoForHintFromEmptySpace(
        const QPointF& scene_pos) const
{
    auto chain_id = rdkit_extensions::toString(m_chain_type) + "1";
    auto monomer =
        rdkit_extensions::makeMonomer(m_res_name, chain_id, 1, false);
    auto monomer_pos = to_mol_xy(scene_pos);
    auto linkage_start = getDefaultDragStartAPModelName();
    // returned monomer is owned by the calling scope
    return HintFragmentMonomerInfo{std::move(monomer), m_monomer_type,
                                   monomer_pos, linkage_start,
                                   NEW_MONOMER_FROM_DRAG};
}

HintFragmentMonomerInfo AbstractDrawMonomerOrMonomericConnectionSceneTool::
    createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
        const AbstractMonomerItem* const monomer_item,
        const std::string& ap_model_name) const
{
    auto [monomer, monomer_type] = get_monomer_and_type(monomer_item);
    // returned monomer is owned by the calling scope
    auto copy_of_monomer = std::make_unique<RDKit::Atom>(*monomer);
    auto monomer_pos = get_coords_for_monomer(monomer);
    auto monomer_idx = static_cast<int>(monomer->getIdx());
    return HintFragmentMonomerInfo{std::move(copy_of_monomer), monomer_type,
                                   monomer_pos, ap_model_name, monomer_idx};
}

HintFragmentMonomerInfo AbstractDrawMonomerOrMonomericConnectionSceneTool::
    createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
        const AbstractMonomerItem* const monomer_item,
        const UnboundMonomericAttachmentPointItem* const ap_item) const
{
    auto ap_model_name = ap_item->getAttachmentPoint().model_name;
    return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
        monomer_item, ap_model_name);
}

HintFragmentMonomerInfo AbstractDrawMonomerOrMonomericConnectionSceneTool::
    createHintFragmentMonomerInfoForHintToDirection(
        const HintFragmentMonomerInfo& start_monomer_info,
        const Direction direction) const
{
    auto pos =
        get_default_coords_for_bound_monomer(start_monomer_info.pos, direction);
    // We'll determine the correct values for chain id and residue number
    // when/if the structure is actually added to MolModel. For the hint,
    // though, just generate something reasonable looking.
    auto chain_id = rdkit_extensions::toString(m_chain_type) + "1";
    auto res_num = 2;
    // returned monomer is owned by the calling scope
    auto monomer =
        rdkit_extensions::makeMonomer(m_res_name, chain_id, res_num, false);
    auto ap_model_name = get_attachment_point_for_new_monomer(
        start_monomer_info.monomer_type, start_monomer_info.ap_model_name,
        m_monomer_type);
    return HintFragmentMonomerInfo{std::move(monomer), m_monomer_type, pos,
                                   ap_model_name, NEW_MONOMER_FROM_DRAG};
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::onLeftButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonDragStart(event);
    m_drag_start_monomer_item = getTopMonomerItemAt(m_mouse_press_scene_pos);
    if (m_drag_start_monomer_item == nullptr) {
        m_drag_start_ap_model_name = getDefaultDragStartAPModelName();
    } else {
        m_drag_start_ap_model_name =
            getAPModelNameAtPosOverMonomer(m_mouse_press_scene_pos);
    }
    // clear the attachment point labels. Note that this must happen *after* the
    // getUnboundAttachmentPointAt call, since that method requires the
    // attachment points to be in the Scene
    clearAttachmentPointsLabelsAndHintFragmentItem();

    auto [drag_end_info, _] = getDragEndInfo(event->scenePos());
    auto handled = createDragHintIfDragStartValid(drag_end_info);
    if (handled) {
        // hide the cursor hint while the hint structure is shown. Note that the
        // cursor hint will have already been hidden if we started this drag
        // from a monomer, since hovering over a monomer hides the cursor hint
        // (since it shows the hint structure).  If this drag started from empty
        // space, though, then this call is necessary to hide the cursor hint.
        setCursorHintShown(false);
    } else {
        m_drag_start_monomer_item = nullptr;
    }
    m_drag_ignored = !handled;
}

std::string AbstractDrawMonomerOrMonomericConnectionSceneTool::
    getAPModelNameAtPosOverMonomer(const QPointF& scene_pos) const
{
    auto ap_item = getUnboundAttachmentPointAt(scene_pos, false);
    if (ap_item == nullptr) {
        // this monomer has no available attachment points
        return "";
    } else {
        return ap_item->getAttachmentPoint().model_name;
    }
}
bool AbstractDrawMonomerOrMonomericConnectionSceneTool::
    createDragHintIfDragStartValid(const DragEndInfo& drag_end_info)
{
    clearHintFragmentItem();
    auto hint_start_monomer_info = getHintFragmentMonomerInfoForDragStart();
    if (!hint_start_monomer_info.has_value()) {
        return false;
    }
    auto hint_end_monomer_info = getHintFragmentMonomerInfoForDragEnd(
        *hint_start_monomer_info, drag_end_info);
    if (!hint_end_monomer_info.has_value()) {
        return false;
    }
    createHintFragmentItem(*hint_start_monomer_info, *hint_end_monomer_info);

    return true;
}

std::pair<DragEndInfo, AbstractMonomerItem*>
AbstractDrawMonomerOrMonomericConnectionSceneTool::getDragEndInfo(
    const QPointF& scene_pos)
{
    auto* hovered_monomer_item = getTopMonomerItemAt(scene_pos);
    if (hovered_monomer_item == m_drag_start_monomer_item) {
        // we can't drag from a monomer to itself
        hovered_monomer_item = nullptr;
    }
    UnboundMonomericAttachmentPointItem* drag_end_ap_item =
        (hovered_monomer_item == nullptr)
            ? nullptr
            : getUnboundDragEndAttachmentPointAt(scene_pos);

    // if we're dragging between two existing monomers, make sure that we'd be
    // able to actually add a connection between them. Any pair of monomers can
    // have at most one standard backbone connection and one custom connection
    // (i.e. not a standard backbone connection) between them, since our RDKit
    // monomer format currently doesn't support more than that
    if (m_drag_start_monomer_item != nullptr && drag_end_ap_item != nullptr) {
        auto ap_name_end = drag_end_ap_item->getAttachmentPoint().model_name;
        if (!dragCanFormConnectionTo(hovered_monomer_item, ap_name_end)) {
            drag_end_ap_item = nullptr;
        };
    }

    DragEndInfo drag_end_info;
    if (drag_end_ap_item == nullptr) {
        // there's no available attachment point at the cursor (or there is one
        // but we can't use it because this pair of monomers already has that
        // type of connection), so we display the drag hint in a direction
        drag_end_info = getDragEndInfoForNonMonomerPos(scene_pos);
    } else {
        drag_end_info = std::make_pair(hovered_monomer_item, drag_end_ap_item);
    }
    return {drag_end_info, hovered_monomer_item};
}

bool AbstractDrawMonomerOrMonomericConnectionSceneTool::dragCanFormConnectionTo(
    const AbstractMonomerItem* const hovered_monomer_item,
    const std::string& ap_name_end) const
{
    auto* monomer_start = m_drag_start_monomer_item->getAtom();
    auto* monomer_end = hovered_monomer_item->getAtom();
    auto* connection = m_mol_model->getMol()->getBondBetweenAtoms(
        monomer_start->getIdx(), monomer_end->getIdx());
    if (connection != nullptr) {
        if (ap_name_end == H_BOND_AP_MODEL_NAME ||
            is_hydrogen_bond(connection)) {
            // If either the new bond or the existing bond is a hydrogen bond,
            // then we're trying to mix a hydrogen bond and a covalent bond in
            // the same connection, which isn't allowed. If both of them are
            // hydrogen bonds, then we're trying to create a hydrogen bond that
            // already exists, which isn't allowed.
            return false;
        }
        auto [linkage, flipped] =
            build_linkage_string(m_drag_start_ap_model_name, ap_name_end);
        bool is_custom_bond =
            get_is_custom_bond(monomer_start, monomer_end, linkage);
        if (is_custom_bond && connection->hasProp(CUSTOM_BOND)) {
            // this would be a custom connection (i.e. not a standard
            // backbone connection) and this monomer pair already has one
            return false;
        } else if (!is_custom_bond &&
                   (!connection->hasProp(CUSTOM_BOND) ||
                    contains_two_monomer_linkages(connection))) {
            // this would be a standard backbone connection and this monomer
            // pair already has one (presumably in the opposite direction,
            // since otherwise the attachment points would have already been
            // bound)
            return false;
        }
    }
    return true;
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonDragMove(event);
    if (m_drag_ignored) {
        return;
    }
    auto [drag_end_info, hovered_monomer_item] =
        getDragEndInfo(event->scenePos());

    if (hovered_monomer_item != m_drag_end_monomer_item) {
        if (m_drag_end_monomer_item) {
            clearDragEndAttachmentPointsLabels();
        }
        m_drag_end_monomer_item = hovered_monomer_item;
        if (hovered_monomer_item != nullptr) {
            labelAttachmentPointsOnDragEndMonomer(
                hovered_monomer_item->getAtom(), hovered_monomer_item);
        }
    }

    if (drag_end_info != m_drag_end_info) {
        m_drag_end_info = drag_end_info;
        // we know the drag start was valid, since otherwise m_drag_ignored
        // would have been true
        createDragHintIfDragStartValid(drag_end_info);

        // update which unbound attachment point is highlighted
        if (std::holds_alternative<MonomerAndAPItems>(drag_end_info)) {
            auto active_ap_item =
                std::get<MonomerAndAPItems>(drag_end_info).second;
            for (auto* ap_item : m_drag_end_unbound_ap_items) {
                auto color = ap_item == active_ap_item
                                 ? STRUCTURE_HINT_COLOR
                                 : m_drag_end_inactive_ap_color;
                ap_item->setColor(color);
            }
        }
    }
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    AbstractMonomerSceneTool::onLeftButtonDragRelease(event);
    if (m_drag_ignored) {
        return;
    }

    // we need to figure out what attachment point we're over before we delete
    // the attachment point graphics items
    auto hint_start_monomer_info = getHintFragmentMonomerInfoForDragStart();
    auto [drag_end_info, hovered_monomer_item] =
        getDragEndInfo(event->scenePos());
    // we know that hint_start_monomer_info can't be std::nullopt, since
    // otherwise m_drag_ignored would be true and we would've returned already
    auto hint_end_monomer_info = getHintFragmentMonomerInfoForDragEnd(
        *hint_start_monomer_info, drag_end_info);

    // delete the attachment point graphics items before we modify the structure
    // so that we don't have to worry about monomer graphics items being deleted
    // and automatically destroying their children
    clearAttachmentPointsLabelsAndHintFragmentItem();
    clearDragEndAttachmentPointsLabels();
    m_drag_start_monomer_item = nullptr;
    m_drag_start_ap_model_name.clear();
    m_drag_end_monomer_item = nullptr;
    m_drag_end_info = std::monostate{};

    // now that everything is cleaned up, we can actually add the monomers and
    // connection to MolModel, assuming the drag end represents a valid
    // structure
    if (hint_end_monomer_info.has_value()) {
        addDragStructureToMolModel(*hint_start_monomer_info,
                                   *hint_end_monomer_info);
    }
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::
    addDragStructureToMolModel(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const HintFragmentMonomerInfo& hint_end_monomer_info)
{
    auto add_monomer_to_mol_model_if_new =
        [this](const HintFragmentMonomerInfo& monomer_info) {
            if (monomer_info.atom_idx >= 0) {
                // the monomer already exists in MolModel
                return static_cast<unsigned int>(monomer_info.atom_idx);
            } else {
                auto monomer_idx = m_mol_model->getMol()->getNumAtoms();
                m_mol_model->addMonomer(m_res_name, m_chain_type,
                                        monomer_info.pos);
                return monomer_idx;
            }
        };

    auto undo_raii = m_mol_model->createUndoMacro("Add monomeric connection");
    auto start_monomer_idx =
        add_monomer_to_mol_model_if_new(hint_start_monomer_info);
    auto end_monomer_idx =
        add_monomer_to_mol_model_if_new(hint_end_monomer_info);
    auto mol = m_mol_model->getMol();
    m_mol_model->addMonomericConnection(mol->getAtomWithIdx(start_monomer_idx),
                                        hint_start_monomer_info.ap_model_name,
                                        mol->getAtomWithIdx(end_monomer_idx),
                                        hint_end_monomer_info.ap_model_name);
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::
    labelAttachmentPointsOnDragEndMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item)
{
    labelUnboundAttachmentPointsOnMonomer(
        monomer, monomer_item, m_drag_end_attachment_point_labels_group,
        m_drag_end_unbound_ap_items);
}

void AbstractDrawMonomerOrMonomericConnectionSceneTool::
    clearDragEndAttachmentPointsLabels()
{
    clear_graphics_item_group(m_drag_end_attachment_point_labels_group);
    delete_all_unbound_ap_items(m_drag_end_unbound_ap_items);
}

} // namespace sketcher
} // namespace schrodinger
