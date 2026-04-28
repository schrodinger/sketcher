#pragma once

#include <tuple>
#include <vector>
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"
#include "schrodinger/sketcher/tool/abstract_monomer_scene_tool.h"

namespace RDKit
{
class Atom;
class Bond;
class ROMol;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class AbstractMonomerItem;
class DrawMonomerFragmentSceneTool;
class UnboundMonomericAttachmentPointItem;

/**
 * A description of one possible way that a monomeric fragment can connect to an
 * existing structure: a connection is allowed between the struc_ap_model_name
 * attachment point of a struc_monomer_type monomer in the existing structure,
 * and the frag_monomer_ap_model_name attachment point of the monomer with index
 * frag_monomer_idx in the fragment.
 */
struct MonomerFragmentAttachmentInfo {
    MonomerType struc_monomer_type;
    std::string struc_ap_model_name;
    unsigned int frag_monomer_idx;
    std::string frag_monomer_ap_model_name;
};

/**
 * @return a scene tool for drawing a nucleotide fragment with
 * the given residue names.
 */
std::shared_ptr<DrawMonomerFragmentSceneTool>
get_nucleotide_fragment_scene_tool(const std::string& sugar,
                                   const std::string& base,
                                   const std::string& phos, const Fonts& fonts,
                                   Scene* scene, MolModel* mol_model);

/**
 * A scene tool for drawing monomeric fragments
 */
class SKETCHER_API DrawMonomerFragmentSceneTool
    : public AbstractMonomerSceneTool
{
  public:
    /**
     * @param mol The fragment to draw
     * @param attachment_info A list of all allowable ways to connect this
     * fragment to the existing structure. If the user hovers over a monomer or
     * attachment point on this list, a hint fragment will be drawn showing the
     * bound fragment. If the user clicks, the fragment and connection will be
     * added to the structure. If the same struc_monomer_type appears more than
     * once on this list, the first attachment info will take precedence when
     * the user hovers over the monomer. (E.g., if this list contains an entry
     * for connecting to the C-terminus of a peptide as well as an entry for
     * connecting to the N-terminus of a peptide, then whichever one is listed
     * first will take effect when the user hovers over the middle of a peptide
     * monomer.)
     * @param index_to_center_on_click When the user clicks on an empty space in
     * the Sketcher workspace, the monomer with this index will be placed at the
     * clicked coordinates.
     * @param drag_angle_adjustment When the user clicks in empty space and
     * drags, the hint fragment will be rotated. Dragging in this direction (as
     * measured in radians clockwise from the positive X direction) will count
     * as no rotation. For example, if this angle is pi / 2 (as it is for
     * nucleotides), then clicking and dragging directly down will have no
     * effect on the rotation. Clicking and dragging in any other direction will
     * rotate the nucleobase (which starts out directly below the sugar) in that
     * direction.
     * @param fonts The fonts to use for rendering
     * @param scene The scene that this tool will be used with
     * @param mol_model The MolModel used in scene
     */
    DrawMonomerFragmentSceneTool(
        const RDKit::ROMol& mol,
        const std::vector<MonomerFragmentAttachmentInfo> attachment_info,
        const int index_to_center_on_click, const double drag_angle_adjustment,
        const Fonts& fonts, Scene* scene, MolModel* mol_model);

    // Reimplemented AbstractSceneTool methods
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonPress(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonRelease(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;

  protected:
    RDKit::ROMol m_frag;
    std::vector<MonomerFragmentAttachmentInfo> m_attachment_info;
    int m_index_to_center_on_click;
    bool m_press_was_over_monomer = false;
    double m_drag_angle_adjustment;

    /**
     * @return the attachment point graphics item that we should attach the
     * fragment to when the user hovers or clicks at the specified coordinates.
     * This attachment point must match one of the attachment_info entries
     * passed into the constructor. If no such attachment point exists, nullptr
     * is returned.
     */
    UnboundMonomericAttachmentPointItem* getConnectableUnboundAttachmentPointAt(
        const QPointF& scene_pos,
        const AbstractMonomerItem* hovered_monomer) const;

    /**
     * Clear any existing fragment hint and draw a hint showing the fragment
     * bound to the specified monomer and attachment point.  Note that, if
     * hovered_ap_item is nullptr, no new hint will be drawn.
     */
    void drawFragmentHintFor(
        const AbstractMonomerItem* hovered_item,
        const UnboundMonomericAttachmentPointItem* hovered_ap_item);

    /**
     * @return A copy of the fragment, positioned correctly to interact with the
     * specified monomer and attachment point via the binding described in
     * attachment_info.
     */
    std::shared_ptr<RDKit::RWMol>
    getPositionedFragment(const RDKit::Atom* const hovered_monomer,
                          const UnboundAttachmentPoint& hovered_ap,
                          const MonomerFragmentAttachmentInfo& attachment_info);

    void addBoundFragmentToMolModel(
        const AbstractMonomerItem* hovered_monomer_item,
        const UnboundMonomericAttachmentPointItem* hovered_ap_item);

    void addUnboundFragmentToMolModel(const QPointF& scene_pos,
                                      const double rotation = 0.0);

    /**
     * @return the angle, in radians, that the drag hint should be rotated when
     * the user has dragged to the specified coordinates.
     */
    double getDragAngle(const QPointF& scene_pos);
};

} // namespace sketcher
} // namespace schrodinger
