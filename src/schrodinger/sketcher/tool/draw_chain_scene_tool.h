#pragma once

#include <rdkit/GraphMol/ROMol.h>

#include <QList>
#include <QGraphicsItem>
#include <QGraphicsItemGroup>
#include <QGraphicsPathItem>
#include <QGraphicsPathItem>
#include <QGraphicsSimpleTextItem>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A blue hint to show where bonds would be added if the user released their
 * mouse drag.  Also displays a text label with the number of bonds that would
 * be added.
 */
class HintChainItem : public QGraphicsItemGroup
{
  public:
    HintChainItem(QGraphicsItem* parent = nullptr);

    /**
     * Update the hint to display bonds between all provided atomic coordinates.
     * @param coords Atomic positions in Scene coordinates.
     */
    void setCoords(QList<QPointF> coords);

  protected:
    QGraphicsPathItem m_bonds_item;
    QGraphicsSimpleTextItem m_label_item;
};

/**
 * A scene tool for drawing a non-branched chain of carbons all at once
 */
class SKETCHER_API DrawChainSceneTool : public StandardSceneToolBase
{
  public:
    DrawChainSceneTool(Scene* scene, MolModel* mol_model);

    // Overridden AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    void onLeftButtonDragStart(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;
    QPixmap createDefaultCursorPixmap() const override;

  protected:
    HintChainItem m_hint_chain_item = HintChainItem();

    /**
     * Return info about where the drag started
     * @return A tuple of
     *   - The coordinates for the start of the drag.  If the mouse was pressed
     *     over an atom, then this will be the coordinates of the atom.  If not,
     *     then this will be the mouse down coordinates.
     *   - The atom, if any, that the mouse was pressed over.  nullptr
     *     otherwise.
     */
    std::pair<QPointF, const RDKit::Atom*> getStartPosAndAtom();
};

/**
 * @return A list of atomic positions in Scene coordinates that is
 * appropriate for a mouse cursor drag from start to end
 */
SKETCHER_API QList<QPointF> get_bond_chain_atom_coords(QPointF start,
                                                       QPointF end);

} // namespace sketcher
} // namespace schrodinger
