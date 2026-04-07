#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"

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

/**
 * Information about one monomer of the hint structure fragment (i.e. the blue
 * structure that shows up when the user hovers over a monomer or
 * click-and-drags from a monomer)
 */
struct HintFragmentMonomerInfo {
    RDKit::Atom* monomer;
    MonomerType monomer_type;
    RDGeom::Point3D pos;
    // the model name of this monomer's attachment point that's connected to the
    // bond (e.g. "R2", not "C")
    std::string ap_model_name;
    // the atom index for this monomer in the MolModel molecule. Will be
    // NEW_MONOMER_FROM_DRAG if the monomer is not present in MolModel.
    int atom_idx;
};

constexpr int NEW_MONOMER_FROM_DRAG = -1;

DrawMonomerSceneTool::DrawMonomerSceneTool(
    const std::string& res_name, const rdkit_extensions::ChainType chain_type,
    const Fonts& fonts, Scene* scene, MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model),
    m_res_name(res_name),
    m_chain_type(chain_type),
    m_fonts(fonts)
{
    // we handle predictive highlighting manually in onMouseMove, so we disable
    // StandardSceneToolsBase's predictive highlighting
    m_highlight_types = InteractiveItemFlag::NONE;
    // make sure that the cursor hint font is more easily readable at small
    // size. (Note that m_fonts is a copy of the Scene's fonts, not a reference,
    // so this change won't affect anything else.)
    m_fonts.m_main_label_font.setBold(true);
    m_fonts.updateFontMetrics();
    m_attachment_point_labels_group.setZValue(static_cast<qreal>(ZOrder::HINT));
    if (chain_type == rdkit_extensions::ChainType::PEPTIDE) {
        m_monomer_type = MonomerType::PEPTIDE;
    } else if (chain_type == rdkit_extensions::ChainType::CHEM) {
        m_monomer_type = MonomerType::CHEM;
    } else {
        m_monomer_type = get_na_monomer_type_from_res_name(res_name);
    }
}

DrawMonomerSceneTool::~DrawMonomerSceneTool()
{
    // explicitly erase any attachment point labels. Without this, the bound
    // attachment point labels would be deleted regardless when
    // m_attachment_point_labels_group and
    // m_drag_end_attachment_point_labels_group are destroyed, but the unbound
    // attachment point labels are parented to their monomer, so they would
    // outlive the scene tool without this call.
    clearAttachmentPointsLabelsAndHintFragmentItem();
    clearDragEndAttachmentPointsLabels();
}

std::vector<QGraphicsItem*> DrawMonomerSceneTool::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
    items.push_back(&m_attachment_point_labels_group);
    items.push_back(&m_drag_end_attachment_point_labels_group);
    return items;
}

void DrawMonomerSceneTool::updateColorsAfterBackgroundColorChange(
    bool is_dark_mode)
{
    StandardSceneToolBase::updateColorsAfterBackgroundColorChange(is_dark_mode);
    m_monomer_background_color =
        is_dark_mode ? DARK_BACKGROUND_COLOR : LIGHT_BACKGROUND_COLOR;
    m_unbound_ap_label_color =
        is_dark_mode ? UNBOUND_AP_LABEL_COLOR_DARK_BG : UNBOUND_AP_LABEL_COLOR;
    m_bound_ap_label_color =
        is_dark_mode ? BOUND_AP_LABEL_COLOR_DARK_BG : BOUND_AP_LABEL_COLOR;
    m_drag_end_inactive_ap_color =
        is_dark_mode ? DRAG_END_INACTIVE_AP_DARK_BG : DRAG_END_INACTIVE_AP;
}

AbstractMonomerItem*
DrawMonomerSceneTool::getTopMonomerItemAt(const QPointF& scene_pos) const
{
    // we're only interested in monomers that are part of the actual Sketcher
    // structure, so we need to ignore any graphics items that belong to our
    // hint structure
    auto is_item_part_of_hint_fragment = [this](QGraphicsItem* item) {
        return m_hint_fragment_item != nullptr &&
               (item == m_hint_fragment_item ||
                item->group() == m_hint_fragment_item);
    };

    // check to see if we're over a monomer, monomeric connector, or unbound
    // attachment point item
    for (auto* item : m_scene->items(scene_pos)) {
        if (is_item_part_of_hint_fragment(item)) {
            // ignore the fragment hint that was drawn by this class
            continue;
        }
        if (item_matches_type_flag(item, InteractiveItemFlag::MONOMER)) {
            return static_cast<AbstractMonomerItem*>(item);
        }
        auto* ap_item =
            qgraphicsitem_cast<UnboundMonomericAttachmentPointItem*>(item);
        if (ap_item && ap_item->withinHoverArea(scene_pos)) {
            return static_cast<AbstractMonomerItem*>(item->parentItem());
        }
    }

    // if we're not over any of those, check to see if we're near a monomer. If
    // we are, check to see whether we'd be over one of its attachment points
    // once they're drawn
    QPainterPath near_scene_pos;
    // We want a circle radius that will definitely include the relevant monomer
    // if we're over an attachment point hover area, but we don't care if it's
    // too large since we'll do further filtering below.  Instead of trying to
    // take the font size, etc. into account, we just pick a number that's close
    // to correct and then double it.
    auto radius =
        2 * std::max(UNBOUND_AP_LINE_LENGTH, UNBOUND_AP_MIN_HOVER_HALF_WIDTH);
    near_scene_pos.addEllipse(scene_pos, radius, radius);
    for (auto* item : m_scene->items(near_scene_pos)) {
        if (!item_matches_type_flag(item, InteractiveItemFlag::MONOMER) ||
            is_item_part_of_hint_fragment(item)) {
            continue;
        }
        auto* monomer_item = static_cast<AbstractMonomerItem*>(item);
        auto local_pos = monomer_item->mapFromScene(scene_pos);
        auto* monomer = monomer_item->getAtom();
        auto [bound_aps, unbound_aps] =
            get_attachment_points_for_monomer(monomer);
        for (auto cur_unbound_ap : unbound_aps) {
            auto unbound_ap_hover_area =
                get_hover_area_for_unbound_monomer_attachment_point_item(
                    cur_unbound_ap, monomer_item, m_fonts);
            if (unbound_ap_hover_area.contains(local_pos)) {
                return monomer_item;
            }
        }
    }
    return nullptr;
}

/**
 * @return the unbound attachment point graphics item representing the
 * attachment point with the "best" number. "Best" is defined using the
 * preferred_order list, with earlier numbers in the list preferred over later
 * numbers. Will return nullptr if no attachment points on the preferred_order
 * list are found.
 */
[[nodiscard]] static UnboundMonomericAttachmentPointItem*
find_preferred_attachment_point_by_num(
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

/**
 * Return the unbound attachment point graphics item that represents the
 * attachment point with the lowest number. An attachment points with a custom
 * name will be returned only if there are no numbered attachment point.
 * @param unbound_ap_items The non-empty list of unbound attachment point
 * graphics item
 */
[[nodiscard]] static UnboundMonomericAttachmentPointItem*
find_min_attachment_point_by_num(
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

/**
 * @return the unbound attachment point graphics item representing an attachment
 * point with the given name. Will return nullptr if no such attachment point is
 * found.
 */
[[nodiscard]] static UnboundMonomericAttachmentPointItem*
find_attachment_point_with_name(
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
                                                   NA_BASE_AP_PAIR);
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
            return NA_BASE_AP_PAIR;
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
DrawMonomerSceneTool::getUnboundAttachmentPointAt(
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

static std::tuple<const RDKit::Atom*, MonomerType>
get_monomer_and_type(const AbstractMonomerItem* const monomer_item)
{
    auto* monomer = monomer_item->getAtom();
    auto monomer_type = get_monomer_type(monomer);
    return {monomer, monomer_type};
}

UnboundMonomericAttachmentPointItem*
DrawMonomerSceneTool::getUnboundDragEndAttachmentPointAt(
    const QPointF& scene_pos) const
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
DrawMonomerSceneTool::getHoveredMonomerAndType() const
{
    if (!item_matches_type_flag(m_hovered_item, InteractiveItemFlag::MONOMER)) {
        throw std::runtime_error("No hovered monomer");
    }
    const auto* monomer_item =
        static_cast<const AbstractMonomerItem*>(m_hovered_item);
    return get_monomer_and_type(monomer_item);
}

UnboundMonomericAttachmentPointItem*
DrawMonomerSceneTool::getDefaultUnboundAttachmentPointForHoveredMonomer(
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

bool DrawMonomerSceneTool::clickShouldMutate(
    const RDKit::Atom* monomer, const MonomerType monomer_type) const
{
    return (monomer_type == m_monomer_type &&
            get_monomer_res_name(monomer) != m_res_name);
}

bool DrawMonomerSceneTool::shouldShowPredictiveHighlighting() const
{

    if (m_hovered_item == nullptr ||
        !item_matches_type_flag(m_hovered_item, InteractiveItemFlag::MONOMER)) {
        return false;
    }
    auto [monomer, monomer_type] = getHoveredMonomerAndType();
    if (clickShouldMutate(monomer, monomer_type)) {
        return true;
    }
    return get_default_attachment_point(monomer_type, m_monomer_type,
                                        m_unbound_ap_items) != nullptr;
}

void DrawMonomerSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag logic is handled in onDragMove
        return;
    }
    QPointF scene_pos = event->scenePos();
    auto* item = getTopMonomerItemAt(scene_pos);
    if (item != m_hovered_item) {
        m_hovered_item = item;
        drawAttachmentPointLabelsFor(item);
        if (shouldShowPredictiveHighlighting()) {
            m_predictive_highlighting_item.highlightItem(item);
        } else {
            m_predictive_highlighting_item.clearHighlightingPath();
        }
    }

    auto* hovered_ap_item = getUnboundAttachmentPointAt(scene_pos, true);
    if (hovered_ap_item != m_hovered_ap_item) {
        // update which attachment point is hovered
        m_hovered_ap_item = hovered_ap_item;
        // hide the hovered "nubbin" and draw a hint fragment instead
        for (auto* ap_item : m_unbound_ap_items) {
            ap_item->setVisible(hovered_ap_item == nullptr ||
                                ap_item != hovered_ap_item);
        }
        drawBoundMonomerHintFor(hovered_ap_item);
    }

    // we want to show either the hint fragment or the cursor hint, but not both
    bool drew_hint_fragment = hovered_ap_item != nullptr;
    if (drew_hint_fragment == m_cursor_hint_shown) {
        emit newCursorHintRequested(
            drew_hint_fragment ? QPixmap() : getDefaultCursorPixmap());
        m_cursor_hint_shown = !drew_hint_fragment;
    }
}

void DrawMonomerSceneTool::drawAttachmentPointLabelsFor(
    QGraphicsItem* const item)
{
    clearAttachmentPointsLabelsAndHintFragmentItem();
    if (item == nullptr) {
        return;
    }
    if (item_matches_type_flag(item, InteractiveItemFlag::MONOMER)) {
        // hovering over a monomer
        auto* monomer_item = static_cast<AbstractMonomerItem*>(item);
        const auto* monomer = monomer_item->getAtom();
        labelAttachmentPointsOnHoveredMonomer(monomer, monomer_item);
    } else if (item_matches_type_flag(item,
                                      InteractiveItemFlag::MONOMER_CONNECTOR)) {
        // hovering over a monomeric connector
        auto* connector_item = qgraphicsitem_cast<MonomerConnectorItem*>(item);
        const auto* connector = connector_item->getBond();
        labelAttachmentPointsOnConnector(
            connector, connector_item->isSecondaryConnection());
    }
}

static RDGeom::Point3D get_coords_for_monomer(const RDKit::Atom* const monomer)
{
    auto& conf = monomer->getOwningMol().getConformer();
    return conf.getAtomPos(monomer->getIdx());
}

/**
 * @return coordinates that are roughly BOND_LENGTH units away from the given
 * monomer in the specified direction. (Monomers are laid out in a grid-like
 * pattern, so diagonal directions will lead to bonds that are slightly longer
 * than BOND_LENGTH.)
 */
static RDGeom::Point3D
get_default_coords_for_bound_monomer(const RDGeom::Point3D monomer_pos,
                                     const Direction dir)
{
    auto offset = rdkit_extensions::direction_to_vector(dir);
    offset *= BOND_LENGTH;
    return monomer_pos + offset;
}

/**
 * @overload takes a monomer instead of monomer coordinates
 */
static RDGeom::Point3D
get_default_coords_for_bound_monomer(const RDKit::Atom* const monomer,
                                     const Direction dir)
{
    auto monomer_pos = get_coords_for_monomer(monomer);
    return get_default_coords_for_bound_monomer(monomer_pos, dir);
}

// should only be called when hovering over a monomer
void DrawMonomerSceneTool::drawBoundMonomerHintFor(
    UnboundMonomericAttachmentPointItem* const ap_item)
{
    clearHintFragmentItem();
    if (ap_item == nullptr) {
        return;
    }
    const auto* monomer_item =
        static_cast<const AbstractMonomerItem*>(m_hovered_item);
    auto hovered_monomer_info =
        createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            monomer_item, ap_item->getAttachmentPoint().model_name);
    auto direction = ap_item->getAttachmentPoint().direction;
    auto new_monomer_info = createHintFragmentMonomerInfoForHintToDirection(
        hovered_monomer_info, direction);
    createHintFragmentItem(hovered_monomer_info, new_monomer_info);
}

void DrawMonomerSceneTool::createHintFragmentItem(
    const HintFragmentMonomerInfo& monomer_one_info,
    const HintFragmentMonomerInfo& monomer_two_info)
{
    m_frag = std::make_shared<RDKit::RWMol>();
    m_frag->setProp(HELM_MODEL, true);

    // create the two monomers
    auto first_idx = m_frag->addAtom(monomer_one_info.monomer, true, true);
    auto second_idx = m_frag->addAtom(monomer_two_info.monomer, true, true);

    // create the connection between them
    auto linkage = fmt::format("{}-{}", monomer_one_info.ap_model_name,
                               monomer_two_info.ap_model_name);
    rdkit_extensions::addConnection(*m_frag, first_idx, second_idx, linkage);
    auto bond_index_to_label =
        m_frag->getBondBetweenAtoms(first_idx, second_idx)->getIdx();

    // flag the atoms as monomeric
    for (auto* atom : m_frag->atoms()) {
        set_atom_monomeric(atom);
    }

    // Add a conformer with the atom coordinates
    auto* frag_conf = new RDKit::Conformer(m_frag->getNumAtoms());
    frag_conf->set3D(false);
    frag_conf->setAtomPos(first_idx, monomer_one_info.pos);
    frag_conf->setAtomPos(second_idx, monomer_two_info.pos);
    m_frag->addConformer(frag_conf, true);

    // hide the monomers that already exist in the Scene
    std::vector<size_t> atom_indices_to_hide;
    if (monomer_one_info.atom_idx >= 0) {
        atom_indices_to_hide.push_back(first_idx);
    }
    if (monomer_two_info.atom_idx >= 0) {
        atom_indices_to_hide.push_back(second_idx);
    }

    m_hint_fragment_item = new MonomerHintFragmentItem(
        m_frag, m_fonts, atom_indices_to_hide, bond_index_to_label,
        m_monomer_background_color);
    m_scene->addItem(m_hint_fragment_item);
}

std::string DrawMonomerSceneTool::getDefaultDragStartAPModelName() const
{
    return m_monomer_type == MonomerType::NA_BASE
               ? NA_BASE_AP_PAIR
               : ap_model_name_for(PeptideAP::C);
}

HintFragmentMonomerInfo
DrawMonomerSceneTool::createHintFragmentMonomerInfoForHintFromEmptySpace(
    const QPointF& scene_pos) const
{
    auto chain_id = rdkit_extensions::toString(m_chain_type) + "1";
    auto monomer =
        rdkit_extensions::makeMonomer(m_res_name, chain_id, 1, false);
    auto monomer_pos = to_mol_xy(scene_pos);
    auto linkage_start = getDefaultDragStartAPModelName();
    // returned monomer is owned by the calling scope
    return HintFragmentMonomerInfo(monomer.release(), m_monomer_type,
                                   monomer_pos, linkage_start,
                                   NEW_MONOMER_FROM_DRAG);
}

HintFragmentMonomerInfo DrawMonomerSceneTool::
    createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
        const AbstractMonomerItem* const monomer_item,
        const std::string& ap_model_name) const
{
    auto [monomer, monomer_type] = get_monomer_and_type(monomer_item);
    // returned monomer is owned by the calling scope
    auto copy_of_monomer = new RDKit::Atom(*monomer);
    auto monomer_pos = get_coords_for_monomer(monomer);
    return HintFragmentMonomerInfo(copy_of_monomer, monomer_type, monomer_pos,
                                   ap_model_name, monomer->getIdx());
}

HintFragmentMonomerInfo DrawMonomerSceneTool::
    createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
        const AbstractMonomerItem* const monomer_item,
        const UnboundMonomericAttachmentPointItem* const ap_item) const
{
    auto ap_model_name = ap_item->getAttachmentPoint().model_name;
    return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
        monomer_item, ap_model_name);
}

HintFragmentMonomerInfo
DrawMonomerSceneTool::createHintFragmentMonomerInfoForHintToDirection(
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
    return HintFragmentMonomerInfo(monomer.release(), m_monomer_type, pos,
                                   ap_model_name, NEW_MONOMER_FROM_DRAG);
}

void DrawMonomerSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonClick(event);
    QPointF scene_pos = event->scenePos();
    auto* item = getTopMonomerItemAt(scene_pos);

    if (item == nullptr) {
        // the click was on empty space, so create a new monomer here
        auto mol_pos = to_mol_xy(scene_pos);
        m_mol_model->addMonomer(m_res_name, m_chain_type, mol_pos);
    } else {
        auto [monomer, monomer_type] = get_monomer_and_type(item);
        std::optional<UnboundAttachmentPoint> clicked_ap;
        auto ap_item = getUnboundAttachmentPointAt(scene_pos, true);
        if (ap_item != nullptr) {
            clicked_ap = ap_item->getAttachmentPoint();
        }

        if (clicked_ap.has_value()) {
            // the user clicked on an attachment point or this monomer has a
            // default attachment point for this tool, so add a new monomer
            // bound to that attachment point
            auto new_monomer_ap_name = get_attachment_point_for_new_monomer(
                monomer_type, clicked_ap->model_name, m_monomer_type);
            auto new_pos = get_default_coords_for_bound_monomer(
                monomer, clicked_ap->direction);
            // the attachment point labels won't be valid once the new monomer
            // is added, so clear them now (otherwise we risk a crash)
            clearAttachmentPointsLabelsAndHintFragmentItem();
            m_mol_model->addBoundMonomer(m_res_name, m_chain_type, new_pos,
                                         new_monomer_ap_name, monomer,
                                         clicked_ap->model_name);

        } else if (clickShouldMutate(monomer, monomer_type)) {
            // the user clicked directly on the monomer and the clicked
            // monomer's residue name is different than the tool's, so we mutate
            // the clicked monomer
            clearAttachmentPointsLabelsAndHintFragmentItem();
            m_mol_model->mutateMonomers({monomer}, m_res_name, m_monomer_type);
        }
    }
}

void DrawMonomerSceneTool::onLeftButtonDragStart(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragStart(event);
    m_drag_start_monomer_item = getTopMonomerItemAt(m_mouse_press_scene_pos);
    if (m_drag_start_monomer_item == nullptr) {
        m_drag_start_ap_model_name = getDefaultDragStartAPModelName();
    } else {
        auto ap_item =
            getUnboundAttachmentPointAt(m_mouse_press_scene_pos, false);
        if (ap_item == nullptr) {
            // this monomer has no available attachment points. In this
            // scenario, createDragHint will return false below and we'll ignore
            // the drag.
            m_drag_start_ap_model_name = "";
        } else {
            m_drag_start_ap_model_name =
                ap_item->getAttachmentPoint().model_name;
        }
    }
    // clear the attachment point labels. Note that this must happen *after* the
    // getUnboundAttachmentPointAt call, since that method requires the
    // attachment points to be in the Scene
    clearAttachmentPointsLabelsAndHintFragmentItem();

    // since the drag just started, we can safely assume that we're not over a
    // different monomer than m_drag_start_monomer_item, which means we want a
    // drag to direction, not to another's monomer attachment point
    auto dir = getDragDirection(event->scenePos());
    auto handled = createDragHintIfDragStartValid(dir);
    if (!handled) {
        m_drag_start_monomer_item = nullptr;
    }
    m_drag_ignored = !handled;
}

bool DrawMonomerSceneTool::createDragHintIfDragStartValid(
    const DragEndInfo& drag_end_info)
{
    clearHintFragmentItem();
    auto hint_start_monomer_info = getHintFragmentMonomerInfoForDragStart();
    if (!hint_start_monomer_info.has_value()) {
        return false;
    }
    HintFragmentMonomerInfo hint_end_monomer_info =
        getHintFragmentMonomerInfoForDragEnd(*hint_start_monomer_info,
                                             drag_end_info);
    createHintFragmentItem(*hint_start_monomer_info, hint_end_monomer_info);

    return true;
}

std::optional<HintFragmentMonomerInfo>
DrawMonomerSceneTool::getHintFragmentMonomerInfoForDragStart()
{
    if (m_drag_start_monomer_item == nullptr) {
        // the drag started over empty space
        if (m_monomer_type == MonomerType::PEPTIDE ||
            m_monomer_type == MonomerType::NA_BASE) {
            return createHintFragmentMonomerInfoForHintFromEmptySpace(
                m_mouse_press_scene_pos);
        } else {
            // it doesn't make much biological sense to connect this monomer
            // type to itself, which is what would typically happen when
            // dragging from empty space, so don't start the drag
            return std::nullopt;
        }
    } else if (!m_drag_start_ap_model_name.empty()) {
        // the drag started over a monomer and it has an available attachment
        // point
        return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            m_drag_start_monomer_item, m_drag_start_ap_model_name);
    } else {
        // the drag started over a monomer, but that monomer has no available
        // unbound attachment points so we can't drag from it
        return std::nullopt;
    }
}

HintFragmentMonomerInfo
DrawMonomerSceneTool::getHintFragmentMonomerInfoForDragEnd(
    const HintFragmentMonomerInfo& hint_start_monomer_info,
    const DragEndInfo& drag_end_info)
{
    if (std::holds_alternative<Direction>(drag_end_info)) {
        auto dir = std::get<Direction>(drag_end_info);
        return createHintFragmentMonomerInfoForHintToDirection(
            hint_start_monomer_info, dir);
    } else {
        auto [hovered_monomer_item, drag_end_ap_item] =
            std::get<MonomerAndAPItems>(drag_end_info);
        return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            hovered_monomer_item, drag_end_ap_item);
    }
}

std::pair<DragEndInfo, AbstractMonomerItem*>
DrawMonomerSceneTool::getDragEndInfo(const QPointF& scene_pos)
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

    DragEndInfo direction_or_attachment_point;
    if (drag_end_ap_item == nullptr) {
        // there's no available attachment point at the cursor, so we drag the
        // drag hint in a direction
        direction_or_attachment_point = getDragDirection(scene_pos);
    } else {
        direction_or_attachment_point =
            std::make_pair(hovered_monomer_item, drag_end_ap_item);
    }
    return {direction_or_attachment_point, hovered_monomer_item};
}

Direction
DrawMonomerSceneTool::getDragDirection(const QPointF& cur_scene_pos) const
{
    const qreal dx = cur_scene_pos.x() - m_mouse_press_scene_pos.x();
    const qreal dy = cur_scene_pos.y() - m_mouse_press_scene_pos.y();
    const qreal abs_dx = std::fabs(dx);
    const qreal abs_dy = std::fabs(dy);
    const qreal sqrt2_plus_1 = std::sqrt(2.0) + 1.0;

    if (abs_dy * sqrt2_plus_1 <= abs_dx) {
        // Within 22.5 degrees of horizontal
        return dx >= 0 ? Direction::E : Direction::W;
    } else if (abs_dy >= abs_dx * sqrt2_plus_1) {
        // Within 22.5 degrees of vertical (Qt +Y is down = South)
        return dy >= 0 ? Direction::S : Direction::N;
    } else {
        // Diagonal
        if (dx >= 0) {
            return dy >= 0 ? Direction::SE : Direction::NE;
        } else {
            return dy >= 0 ? Direction::SW : Direction::NW;
        }
    }
}

void DrawMonomerSceneTool::onLeftButtonDragMove(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragMove(event);
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

void DrawMonomerSceneTool::onLeftButtonDragRelease(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonDragRelease(event);
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
    // connection to MolModel
    addDragStructureToMolModel(*hint_start_monomer_info, hint_end_monomer_info);
}

void DrawMonomerSceneTool::addDragStructureToMolModel(
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

QPixmap DrawMonomerSceneTool::createDefaultCursorPixmap() const
{
    // the specific number used here (the "1") doesn't matter - we just need any
    // number to form a proper chain ID
    auto chain_id = rdkit_extensions::toString(m_chain_type) + "1";
    auto monomer =
        rdkit_extensions::makeMonomer(m_res_name, chain_id, 1, false);

    std::shared_ptr<AbstractMonomerItem> monomer_item;
    monomer_item.reset(get_monomer_graphics_item(monomer.get(), m_fonts));
    monomer_item->setMonomerColors(Qt::GlobalColor::transparent,
                                   CURSOR_HINT_COLOR, CURSOR_HINT_COLOR);
    // make sure that the cursor hint is at least a little smaller than an
    // actual monomer
    auto min_scene_size =
        CURSOR_HINT_IMAGE_SIZE * MONOMER_CURSOR_HINT_MIN_SCENE_SIZE_SCALE;
    return cursor_hint_from_graphics_item(monomer_item.get(), min_scene_size);
}

void DrawMonomerSceneTool::labelAttachmentPointsOnHoveredMonomer(
    const RDKit::Atom* const monomer, AbstractMonomerItem* const monomer_item)
{
    labelAttachmentPointsOnMonomer(monomer, monomer_item,
                                   m_attachment_point_labels_group,
                                   m_unbound_ap_items);
}
void DrawMonomerSceneTool::labelAttachmentPointsOnDragEndMonomer(
    const RDKit::Atom* const monomer, AbstractMonomerItem* const monomer_item)
{
    labelAttachmentPointsOnMonomer(monomer, monomer_item,
                                   m_drag_end_attachment_point_labels_group,
                                   m_drag_end_unbound_ap_items);
}

void DrawMonomerSceneTool::labelAttachmentPointsOnMonomer(
    const RDKit::Atom* const monomer, AbstractMonomerItem* const monomer_item,
    QGraphicsItemGroup& attachment_point_labels_group,
    std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items)
{
    auto [bound_aps, unbound_aps] = get_attachment_points_for_monomer(monomer);
    for (auto& cur_ap : bound_aps) {
        auto* item = create_label_for_bound_attachment_point(
            monomer, cur_ap.bound_monomer, cur_ap.is_secondary_connection,
            cur_ap.display_name, m_bound_ap_label_color, m_fonts, m_scene);
        if (item != nullptr) {
            attachment_point_labels_group.addToGroup(item);
        }
    }
    for (auto& cur_ap : unbound_aps) {
        auto* item = new UnboundMonomericAttachmentPointItem(
            cur_ap, monomer_item, m_unbound_ap_label_color, m_fonts);
        unbound_ap_items.push_back(item);
    }
}

void DrawMonomerSceneTool::labelAttachmentPointsOnConnector(
    const RDKit::Bond* const connector, const bool is_secondary_connection)
{
    auto ap_label_items = create_attachment_point_labels_for_connector(
        connector, is_secondary_connection, m_bound_ap_label_color, m_fonts,
        m_scene);
    for (auto* item : ap_label_items) {
        m_attachment_point_labels_group.addToGroup(item);
    }
}

template <typename T>
static void clear_graphics_item_group_and_list(QGraphicsItemGroup& group,
                                               std::vector<T*>& items_list)
{
    for (auto* item : group.childItems()) {
        group.removeFromGroup(item);
        delete item;
    }
    for (auto* item : items_list) {
        delete item;
    }
    items_list.clear();
}

void DrawMonomerSceneTool::clearAttachmentPointsLabelsAndHintFragmentItem()
{
    clear_graphics_item_group_and_list(m_attachment_point_labels_group,
                                       m_unbound_ap_items);
    m_hovered_ap_item = nullptr;
    clearHintFragmentItem();
}

void DrawMonomerSceneTool::clearHintFragmentItem()
{
    delete m_hint_fragment_item;
    m_hint_fragment_item = nullptr;
    m_frag.reset();
}

void DrawMonomerSceneTool::clearDragEndAttachmentPointsLabels()
{
    clear_graphics_item_group_and_list(m_drag_end_attachment_point_labels_group,
                                       m_drag_end_unbound_ap_items);
}

} // namespace sketcher
} // namespace schrodinger
