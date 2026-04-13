#pragma once

#include <string>
#include <tuple>
#include <utility>
#include <variant>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

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
struct HintFragmentMonomerInfo;

using MonomerAndAPItems =
    std::pair<AbstractMonomerItem*, UnboundMonomericAttachmentPointItem*>;
using DragEndInfo = std::variant<std::monostate, rdkit_extensions::Direction,
                                 MonomerAndAPItems>;

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
 * A scene tools that draws a monomer
 */
class SKETCHER_API DrawMonomerSceneTool : public StandardSceneToolBase
{
  public:
    DrawMonomerSceneTool(const std::string& res_name,
                         const rdkit_extensions::ChainType chain_type,
                         const Fonts& fonts, Scene* scene, MolModel* mol_model);
    virtual ~DrawMonomerSceneTool();

    // Reimplemented AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragStart(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;
    void updateColorsAfterBackgroundColorChange(bool is_dark_mode) override;

  protected:
    std::string m_res_name;
    rdkit_extensions::ChainType m_chain_type;
    Fonts m_fonts;
    MonomerType m_monomer_type;
    QGraphicsItemGroup m_attachment_point_labels_group;
    const QGraphicsItem* m_hovered_item = nullptr;
    const UnboundMonomericAttachmentPointItem* m_hovered_ap_item = nullptr;
    std::vector<UnboundMonomericAttachmentPointItem*> m_unbound_ap_items;
    MonomerHintFragmentItem* m_hint_fragment_item = nullptr;
    std::shared_ptr<RDKit::RWMol> m_frag = nullptr;
    QColor m_monomer_background_color = LIGHT_BACKGROUND_COLOR;
    QColor m_unbound_ap_label_color = UNBOUND_AP_LABEL_COLOR;
    QColor m_bound_ap_label_color = BOUND_AP_LABEL_COLOR;
    bool m_cursor_hint_shown = true;

    bool m_drag_ignored;
    AbstractMonomerItem* m_drag_start_monomer_item;
    std::string m_drag_start_ap_model_name;
    AbstractMonomerItem* m_drag_end_monomer_item;
    DragEndInfo m_drag_end_info;
    QGraphicsItemGroup m_drag_end_attachment_point_labels_group;
    std::vector<UnboundMonomericAttachmentPointItem*>
        m_drag_end_unbound_ap_items;
    QColor m_drag_end_inactive_ap_color;

    QPixmap createDefaultCursorPixmap() const override;

    /**
     * Label all attachment points on the given hovered monomer
     */
    void labelAttachmentPointsOnHoveredMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item);

    /**
     * Label all attachment points on the given monomer at the end of the
     * current click-and-drag operation
     */
    void labelAttachmentPointsOnDragEndMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item);

    /**
     * Label all attachment points on the given monomer, placing all newly
     * created graphics items into the specified group and list. (You probably
     * want to call either labelAttachmentPointsOnHoveredMonomer or
     * labelAttachmentPointsOnDragEndMonomer instead of this.)
     */
    void labelAttachmentPointsOnMonomer(
        const RDKit::Atom* const monomer,
        AbstractMonomerItem* const monomer_item,
        QGraphicsItemGroup& attachment_point_labels_group,
        std::vector<UnboundMonomericAttachmentPointItem*>& unbound_ap_items);

    /**
     * Label both attachment points for the given monomeric connector
     * @param connector The monomeric connector to label
     * @param is_secondary_connection Whether this method should label the
     * secondary connection of the bond instead of the primary connection.
     * Secondary connections occur when a single RDKit::Bond* represents
     * multiple connections, e.g. two neighboring cysteines that are disulfide
     * bonded to each other.
     */
    void labelAttachmentPointsOnConnector(const RDKit::Bond* const connector,
                                          const bool is_secondary_connection);

    /**
     * Label the specified attachment point
     * @param monomer The monomer containing the attachment point to label
     * @param bound_monomer The other monomer involved in the connection
     * @param is_secondary_connection Whether we should label the secondary
     * connection of the bond instead of the primary connection.
     * @param label The text to display in the attachment point label
     */
    void labelBoundAttachmentPoint(const RDKit::Atom* const monomer,
                                   const RDKit::Atom* const bound_monomer,
                                   const bool is_secondary_connection,
                                   const std::string& label);

    void addAttachmentPointLabel(const QString& label,
                                 const QRectF& label_rect);

    /**
     * Clear all attachment point labels drawn on the hovered monomer. Also
     * clear the hint fragment item and the associated structure.
     */
    void clearAttachmentPointsLabelsAndHintFragmentItem();

    /**
     * Remove the hint fragment item from the scene, then destroy it as well as
     * the RDKit structure object it represents.
     */
    void clearHintFragmentItem();

    /**
     * Clear all attachment point labels drawn on the existing monomer at the
     * end of the current click-and-drag operation. Note that this method *does
     * not* clear the hint fragment item.
     */
    void clearDragEndAttachmentPointsLabels();

    /**
     * Return the top graphics item representing a monomer at the given
     * coordinates. If there's no such graphics item, return nullptr.  Note that
     * this method will consider the cursor to be over a monomer if the cursor
     * is over an unbound attachment point belonging to that monomer, or if the
     * cursor *would be* over an unbound attachment point once it's drawn.
     * @param scene_pos The position in Scene coordinates
     */
    AbstractMonomerItem* getTopMonomerItemAt(const QPointF& scene_pos) const;

    /**
     * Clear any existing attachment point labels and draw new ones for the
     * specified monomeric graphics item. If item is nullptr, then all labels
     * will be cleared and no new ones will be drawn.
     */
    void drawAttachmentPointLabelsFor(QGraphicsItem* const item);

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
    UnboundMonomericAttachmentPointItem*
    getDefaultUnboundAttachmentPointForHoveredMonomer(
        const bool no_default_if_click_should_mutate) const;

    /**
     * Draw a hint structure showing a monomer bound to the specified attachment
     * point
     */
    void
    drawBoundMonomerHintFor(UnboundMonomericAttachmentPointItem* const ap_item);

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
    bool clickShouldMutate(const RDKit::Atom* monomer,
                           const MonomerType monomer_type) const;

    /**
     * Create a hint fragment containing the two specified monomers and the
     * connection between them, then add this hint fragment to the scene.
     */
    void createHintFragmentItem(const HintFragmentMonomerInfo& monomer_one,
                                const HintFragmentMonomerInfo& monomer_two);

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
    std::pair<DragEndInfo, AbstractMonomerItem*>
    getDragEndInfo(const QPointF& scene_pos);

    void addDragStructureToMolModel(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const HintFragmentMonomerInfo& hint_end_monomer_info);

    /**
     * If the user has began a "valid" click-and-drag, create and return the
     * HintFragmentMonomerInfo object describing the starting location/monomer
     * of the drag. Otherwise, return std::nullopt. See the
     * createDragHintIfDragStartValid docstring for an explanation of valid
     * versus invalid click-and-drags.
     */
    std::optional<HintFragmentMonomerInfo>
    getHintFragmentMonomerInfoForDragStart();

    /**
     * Create and return the HintFragmentMonomerInfo object describing the end
     * of the current click-and-drag operation.
     */
    HintFragmentMonomerInfo getHintFragmentMonomerInfoForDragEnd(
        const HintFragmentMonomerInfo& hint_start_monomer_info,
        const DragEndInfo& drag_end_info);

    /**
     * @return the direction of the specified position relative to the start of
     * the current click-and-drag operation
     */
    rdkit_extensions::Direction
    getDragDirection(const QPointF& cur_scene_pos) const;
};

} // namespace sketcher
} // namespace schrodinger
