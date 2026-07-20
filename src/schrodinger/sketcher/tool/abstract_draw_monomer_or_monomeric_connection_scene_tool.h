#pragma once

#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <variant>
#include <vector>

#include <rdkit/Geometry/point.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/abstract_monomer_scene_tool.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{

namespace rdkit_extensions
{
enum class ChainType;
enum class Direction;
} // namespace rdkit_extensions

namespace sketcher
{

class MonomerHintFragmentItem;
class UnboundMonomericAttachmentPointItem;
class MonomerConnectorItem;

using MonomerAndAPItems =
    std::pair<AbstractMonomerItem*, UnboundMonomericAttachmentPointItem*>;
using DragEndInfo = std::variant<std::monostate, rdkit_extensions::Direction,
                                 QPointF, MonomerAndAPItems>;

/**
 * Information about one monomer of the hint structure fragment (i.e. the blue
 * structure that shows up when the user hovers over a monomer or
 * click-and-drags from a monomer)
 */
struct HintFragmentMonomerInfo {
    std::unique_ptr<RDKit::Atom> monomer;
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
constexpr int NEW_MONOMER_FROM_INVALID_DRAG = -2;

/**
 * Return the default unbound attachment point; that is, the attachment point
 * that should be selected when the user hovers over a monomer.
 * @param hovered_type The type of monomer being hovered over
 * @param tool_type  The type of monomer that would be drawn by the active scene
 * tool
 * @param unbound_ap_items A list of all graphics items representing unbound
 * attachment points of the hovered monomer
 */
SKETCHER_API UnboundMonomericAttachmentPointItem* get_default_attachment_point(
    const MonomerType hovered_type, const MonomerType tool_type,
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items);

/**
 * When adding a new monomer bound to an existing monomer, determine the
 * appropriate attach point to use for the new-monomer-end of the connection.
 */
SKETCHER_API std::string
get_attachment_point_for_new_monomer(const MonomerType existing_monomer_type,
                                     const std::string_view existing_monomer_ap,
                                     const MonomerType new_monomer_type);

/**
 * @return the monomer represented by the given graphics item, along with its
 * monomer type
 */
SKETCHER_API std::tuple<const RDKit::Atom*, MonomerType>
get_monomer_and_type(const AbstractMonomerItem* const monomer_item);

/**
 * @overload takes a monomer instead of monomer coordinates
 */
SKETCHER_API RDGeom::Point3D
get_default_coords_for_bound_monomer(const RDKit::Atom* const monomer,
                                     const rdkit_extensions::Direction dir);

/**
 * @return the unbound attachment point graphics item representing the
 * attachment point with the "best" number. "Best" is defined using the
 * preferred_order list, with earlier numbers in the list preferred over later
 * numbers. Will return nullptr if no attachment points on the preferred_order
 * list are found.
 */
[[nodiscard]] SKETCHER_API UnboundMonomericAttachmentPointItem*
find_preferred_attachment_point_by_num(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items,
    const std::vector<int>& preferred_order);

/**
 * Return the unbound attachment point graphics item that represents the
 * attachment point with the lowest number. An attachment points with a custom
 * name will be returned only if there are no numbered attachment point.
 * @param unbound_ap_items The non-empty list of unbound attachment point
 * graphics item
 */
[[nodiscard]] SKETCHER_API UnboundMonomericAttachmentPointItem*
find_min_attachment_point_by_num(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items);

/**
 * @return the unbound attachment point graphics item representing an attachment
 * point with the given name. Will return nullptr if no such attachment point is
 * found.
 */
[[nodiscard]] SKETCHER_API UnboundMonomericAttachmentPointItem*
find_attachment_point_with_name(
    const std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items,
    const std::string_view name);

/**
 * An abstract base class for scene tools that draw a monomer or a monomeric
 * connection. This class implements the common behavior for hovering over
 * monomers, drawing hint fragments, and click-and-drag operations. Subclasses
 * must implement getHintFragmentMonomerInfoForDragStart,
 * getHintFragmentMonomerInfoForDragEnd, and clickShouldMutate to specify the
 * tool-specific behavior.
 */
class SKETCHER_API AbstractDrawMonomerOrMonomericConnectionSceneTool
    : public AbstractMonomerSceneTool
{
  public:
    AbstractDrawMonomerOrMonomericConnectionSceneTool(
        const std::string& res_name,
        const rdkit_extensions::ChainType chain_type, const Fonts& fonts,
        Scene* scene, MolModel* mol_model);
    virtual ~AbstractDrawMonomerOrMonomericConnectionSceneTool();

    // Reimplemented AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    void onStructureUpdated() override;
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragStart(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;
    void updateColorsAfterBackgroundColorChange(bool is_dark_mode) override;

  protected:
    std::string m_res_name;
    rdkit_extensions::ChainType m_chain_type =
        rdkit_extensions::ChainType::CHEM;
    MonomerType m_monomer_type;
    Fonts m_bolded_fonts;

    bool m_drag_ignored;
    AbstractMonomerItem* m_drag_start_monomer_item = nullptr;
    std::string m_drag_start_ap_model_name;
    AbstractMonomerItem* m_drag_end_monomer_item = nullptr;
    DragEndInfo m_drag_end_info;
    QGraphicsItemGroup m_drag_end_attachment_point_labels_group;
    std::vector<UnboundMonomericAttachmentPointItem*>
        m_drag_end_unbound_ap_items;
    QColor m_drag_end_inactive_ap_color;

    /**
     * Label all attachment points on the given monomer at the end of the
     * current click-and-drag operation
     */
    void labelAttachmentPointsOnDragEndMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item);

    /**
     * Clear all attachment point labels drawn on the existing monomer at the
     * end of the current click-and-drag operation. Note that this method *does
     * not* clear the hint fragment item.
     */
    void clearDragEndAttachmentPointsLabels();

    /**
     * Return the unbound attachment point graphics item that should be "active"
     * (i.e. that we'd draw a connection for if the user clicked) for the given
     * coordinates. If the coordinates are over an unbound attachment point
     * item, that item will be returned. If the coordinates are over the
     * monomer, then the default unbound attachment point will be returned,
     * assuming one exists for the current tool. If no default attachment point
     * exists (e.g. if the tool and the hovered monomer are of different
     * molecule types), then nullptr will be returned.
     * @param scene_pos the coordinates to check
     * @param no_default_if_click_should_mutate If true, this function will
     * return nullptr if scene_pos is over the monomer itself and clicking on
     * the monomer should trigger a mutation (e.g. if the user has the alanine
     * tool selected and is hovering over a proline).  If false, this function
     * will return the default unbound attachment point in this scenario,
     * assuming one exists. (This value should normally be true if the user is
     * hovering over or has clicked on the monomer, and false if the user has
     * started a drag from the monomer.)
     *
     * @note This method only returns accurate results for the currently hovered
     * monomer, and assumes that the unbound attachment point graphics items
     * have already been drawn for this monomer.
     */
    UnboundMonomericAttachmentPointItem* getUnboundAttachmentPointAt(
        const QPointF& scene_pos,
        const bool no_default_if_click_should_mutate) const;

    /**
     * Return the unbound attachment point graphics item that should be "active"
     * (i.e. that we'd draw a connection for if the user released the mouse
     * button) for the given coordinates, assuming that the user is in the
     * middle of a click-and-drag.  If the coordinates are over the monomer,
     * then the default unbound attachment point will be returned based on the
     * monomer and attachment point that the drag was started from, assuming one
     * exists.
     */
    UnboundMonomericAttachmentPointItem*
    getUnboundDragEndAttachmentPointAt(const QPointF& scene_pos) const;

    /**
     * @return the unbound attachment point that should be active when the user
     * is hovering over the monomer itself
     * @param no_default_if_click_should_mutate If true, this function will
     * return nullptr if clicking on the monomer should trigger a mutation (e.g.
     * if the user has the alanine tool selected and is hovering over a
     * proline).  If false, this function will return the default unbound
     * attachment point in this scenario, assuming one exists. (This value
     * should normally be true if the user is hovering over or has clicked on
     * the monomer, and false if the user has started a drag from the monomer.)
     */
    virtual UnboundMonomericAttachmentPointItem*
    getDefaultUnboundAttachmentPointForHoveredMonomer(
        const bool no_default_if_click_should_mutate) const;

    std::tuple<const RDKit::Atom*, MonomerType>
    getHoveredMonomerAndType() const;

    /**
     * @return whether this tool should show the predictive highlighting outline
     * for the currently hovered monomer.  The highlighting is shown if the
     * hovered monomer can be mutated (i.e. its the same monomer type as this
     * tool but a different residue name) or if there is a reasonable way to
     * connect this tool's monomer to the hovered monomer (an unreasonable
     * connection would be, e.g., connecting a nucleic acid base directly to a
     * phosphate with no intervening sugar.)
     */
    bool shouldShowPredictiveHighlighting() const;

    /**
     * Whether clicking on the specified monomer should mutate that monomer. We
     * only mutate monomers that are the same monomer type as this tool but have
     * a different residue name.
     */
    virtual bool clickShouldMutate(const RDKit::Atom* monomer,
                                   const MonomerType monomer_type) const = 0;

    /**
     * Create a hint fragment containing the two specified monomers and the
     * connection between them, then add this hint fragment to the scene. Note
     * that this function will transfer ownership of the monomers themselves to
     * RDKit.
     */
    virtual void createHintFragmentItem(HintFragmentMonomerInfo& monomer_one,
                                        HintFragmentMonomerInfo& monomer_two);

    /**
     * If the user has started a "valid" click-and-drag operation, create the
     * drag hint and return true. Otherwise, return false. The click-and-drag
     * will be valid unless one of the following happened:
     *  - the user tried to start a drag from an existing monomer that has no
     *    available attachment points
     *  - the user tried to start a drag from empty space while using a nucleic
     *    acid phosphate or sugar tool (It doesn't make any biological sense to
     *    create a dimer of those, so we disallow drags from empty space.)
     * @return
     */
    bool createDragHintIfDragStartValid(const DragEndInfo& drag_end_info);

    /**
     * @return the attachment point to use if user starts a drag from empty
     * space. We only allow drags from empty space for peptides and nucleic acid
     * bases, so this function is only valid for those monomer types. (A drag
     * from empty space normally creates two monomers of the same type, and it
     * doesn't make biological sense to have a dimer of nucleic acid sugars or
     * phosphates.)
     */
    std::string getDefaultDragStartAPModelName() const;

    /**
     * Create a HintFragmentMonomerInfo object describing a new monomer at the
     * given coordinates. Note that the returned monomer is owned by the calling
     * scope.
     */
    HintFragmentMonomerInfo createHintFragmentMonomerInfoForHintFromEmptySpace(
        const QPointF& scene_pos) const;

    /**
     * Create a HintFragmentMonomerInfo object describing the existing monomer
     * and attachment point. Note that the returned monomer is owned by the
     * calling scope.
     */
    HintFragmentMonomerInfo
    createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
        const AbstractMonomerItem* const monomer_item,
        const std::string& ap_model_name) const;

    /**
     * @overload takes a graphics item for the attachment point in place of the
     * model name
     */
    HintFragmentMonomerInfo
    createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
        const AbstractMonomerItem* const monomer_item,
        const UnboundMonomericAttachmentPointItem* const ap_item) const;

    /**
     * Create a HintFragmentMonomerInfo object describing a new monomer at the
     * given coordinates. Note that the returned monomer is owned by the calling
     * scope.
     */
    HintFragmentMonomerInfo createHintFragmentMonomerInfoForHintToDirection(
        const HintFragmentMonomerInfo& start_monomer_info,
        const rdkit_extensions::Direction direction) const;

    /**
     * Determine whether the end of current click-and-drag operation is over an
     * existing monomer or not. If not, return the direction of the drag.
     * @param scene_pos The scene coordinates representing the end of the
     * click-and-drag.
     * @return If scene_pos is over an existing monomer, returns a pair of
     *   - a pair of pointers to the graphics items representing the monomer and
     *     the relevant attachment point. The second value will be nullptr if
     * the monomer does not have any available unbound attachment points.
     *   - a pointer to the monomer's graphics item (i.e. the first pointer of
     *     the pair)
     * If scene_pos is not over an existing monomer, returns a pair of
     *   - the `Direction` enum value representing the direction of the
     *     drag
     *   - nullptr
     */
    virtual std::pair<DragEndInfo, AbstractMonomerItem*>
    getDragEndInfo(const QPointF& scene_pos);

    /**
     * @return Whether the current drag can add a connection to the specified
     * monomer and attachment point model name.
     */
    bool dragCanFormConnectionTo(
        const AbstractMonomerItem* const hovered_monomer_item,
        const std::string& ap_name_end) const;

    virtual void addDragStructureToMolModel(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const HintFragmentMonomerInfo& hint_end_monomer_info);

    /**
     * If the user has began a "valid" click-and-drag, create and return the
     * HintFragmentMonomerInfo object describing the starting location/monomer
     * of the drag. Otherwise, return std::nullopt. See the
     * createDragHintIfDragStartValid docstring for an explanation of valid
     * versus invalid click-and-drags.
     */
    virtual std::optional<HintFragmentMonomerInfo>
    getHintFragmentMonomerInfoForDragStart() = 0;

    /**
     * Create and return the HintFragmentMonomerInfo object describing the end
     * of the current click-and-drag operation, or std::nullopt if the
     * click-and-drag cannot be completed.
     */
    virtual std::optional<HintFragmentMonomerInfo>
    getHintFragmentMonomerInfoForDragEnd(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const DragEndInfo& drag_end_info) = 0;

    /**
     * @return the direction of the specified position relative to the start of
     * the current click-and-drag operation
     */
    rdkit_extensions::Direction
    getDragDirection(const QPointF& cur_scene_pos) const;

    /**
     * Update the graphics items representing unbound attachment points for a
     * monomer when the user hovers over the monomer or its attachment points.
     * @param hovered_ap_item the attachment point that the user is currently
     * hovering over, or the default attachment point if the user is hovering
     * over the monomer itself. May be nullptr if there are no unbound
     * attachment points.
     */
    virtual void updateHoveredUnboundAP(
        UnboundMonomericAttachmentPointItem* hovered_ap_item) = 0;

    /**
     * @return the DragEndInfo for the specified position, which does *not*
     * correspond to a monomer or an attachment point.
     */
    virtual DragEndInfo
    getDragEndInfoForNonMonomerPos(const QPointF& scene_pos) const = 0;

    /**
     * @return the model name of the attachment point to use when the mouse
     * cursor is at the specified coordinates. Returns empty string if there are
     * no unbound attachment points at the coordinates.
     * @param scene_pos The mouse position in scene coordinates, which are
     * guaranteed to be over a monomer or an unbound monomeric attachment point.
     * However, if the coordinates are over a monomer, there's no guarantee that
     * the monomer has any unbound attachment points.
     */
    virtual std::string
    getAPModelNameAtPosOverMonomer(const QPointF& scene_pos) const;
};

} // namespace sketcher
} // namespace schrodinger
