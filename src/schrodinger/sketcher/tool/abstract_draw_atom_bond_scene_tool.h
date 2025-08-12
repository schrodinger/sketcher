#pragma once

#include <tuple>
#include <utility>

#include <QGraphicsLineItem>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * A blue line used to show where a bond would be added if the user clicked or
 * released their mouse drag.
 */
class HintBondItem : public QGraphicsLineItem
{
  public:
    HintBondItem(QGraphicsItem* parent = nullptr);
};

/**
 * A base class for scene tools that draw atoms or bonds (including query atoms
 * and query bonds).
 */
class SKETCHER_API AbstractDrawSceneTool : public StandardSceneToolBase
{
  public:
    AbstractDrawSceneTool(Scene* scene, MolModel* mol_model);

    // Reimplemented AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonRelease(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragStart(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;

  protected:
    HintBondItem m_hint_bond_item = HintBondItem();

    /**
     * Update whether the hint bond is visible
     */
    virtual void setHintBondVisible(bool visible);

    /**
     * Update where the hint bond is drawn.  Note that this will not affect hint
     * bond visibility.
     */
    virtual void updateHintBondPath(const QPointF& start, const QPointF& end);

    /**
     * Modify the specified bond.  Single bonds become double bonds, double
     * bonds become triple bonds, and all other bonds become single bonds.
     */
    void cycleBond(const RDKit::Bond* const bond);

    /**
     * Figure out where we would draw a bond if the user clicked on the
     * specified atom.
     * @param atom The atom to add a bond to
     * @return A tuple of
     *   - The location for the new atom to add (i.e. the atom on the other side
     *     of the bond), given in MolModel coordinates.
     *   - The end point for drawing the hint bond, given in Scene coordinates.
     *   - The existing atom that the default bond should connect to, or a
     *     nullptr if the default bond will require adding a new atom.
     */
    virtual std::tuple<RDGeom::Point3D, QPointF, const RDKit::Atom*>
    getDefaultBondPosition(const RDKit::Atom* const atom) const;

    /**
     * Figure out where we would draw a bond if the user clicked on the
     * specified atom, ignoring any other existing atoms.
     * @param atom The atom to add a bond to
     * @return A tuple of
     *   - The location for the new atom to add (i.e. the atom on the other side
     *     of the bond), given in MolModel coordinates.
     *   - The end point for drawing the hint bond, given in Scene coordinates.
     */
    std::pair<RDGeom::Point3D, QPointF>
    getInitialDefaultBondPosition(const RDKit::Atom* const atom) const;

    /**
     * Determine whether we should start a mouse drag action based on the
     * graphics item (or lack thereof) where the drag started.
     * @return A tuple of
     *   - Whether we should start a mouse drag.
     *   - Where the drag should start from, given in Scene coordinates.
     *   - The atom that the drag was started over, if any.  Nullptr otherwise.
     */
    virtual std::tuple<bool, QPointF, const RDKit::Atom*>
    getDragStartInfo() const;

    /**
     * Determine where we should draw the end of the hint bond during a drag.
     * If the user is currently mousing over the start atom, we'll use the
     * default bond location (i.e. where we would've drawn the bond if the user
     * clicked instead of dragged).  If the user is currently mousing over any
     * other atom, this will be that atom's coordinates.  Otherwise, a bond of
     * the appropriate length is drawn at the nearest 30 degree interval.
     * @param start_pos The start of the drag in Scene coordinates
     * @param start_atom The atom that the cursor was over when starting the
     * drag.  Should be nullptr if the drag didn't start over an atom.
     * @param event_pos The current mouse location in Scene coordinates
     * @return A tuple of
     *   - The location for the end of the hint bond in Scene coordinates.
     *   - The moused-over atom, if any.  Nullptr otherwise.
     */
    virtual std::pair<QPointF, const RDKit::Atom*>
    getBondEndInMousedDirection(const QPointF& start_pos,
                                const RDKit::Atom* const start_atom,
                                const QPointF& event_pos) const;

    /**
     * Determine where we should draw the end of the hint bond during a drag,
     * ignoring any existing atoms.  A bond of the appropriate length will be
     * drawn at the nearest 30 degree interval.
     * @param start The start of the drag in Scene coordinates
     * @param event_pos The current mouse location in Scene coordinates
     * @return The location for the end of the hint bond in Scene coordinates.
     */
    QPointF
    getDefaultBondOffsetInMousedDirection(const QPointF& start,
                                          const QPointF& mouse_pos) const;

    /**
     * Respond to the user clicking on an atom
     * @param atom The atom that was clicked on
     */
    virtual void onAtomClicked(const RDKit::Atom* const atom);

    /**
     * Mutate the specified atom.  Note that the default implementation of this
     * is a no-op.  Subclasses capable of mutating atoms should reimplement this
     * method.
     */
    virtual void mutateAtom(const RDKit::Atom* const atom);

    /**
     * Determine whether we should add a bond if the user were to click on the
     * specified atom.  For example, clicking on a carbon with the "Draw Carbon"
     * tool should add a new bond.  Clicking on a nitrogen with the "Draw
     * Carbon" tool should instead change the nitrogen to a carbon.
     */
    virtual bool
    shouldDrawBondForClickOnAtom(const RDKit::Atom* const atom) const = 0;

    /**
     * Respond to the user clicking on a bond
     * @param bond The bond that was clicked on
     */
    virtual void onBondClicked(const RDKit::Bond* const bond) = 0;

    /**
     * Respond to the user clicking on empty space (no atom or bond)
     * @param pos The location that was clicked on, given in MolModel
     * coordinates
     */
    virtual void onEmptySpaceClicked(const RDGeom::Point3D& pos) = 0;

    /**
     * Add an atom at the specified coordinates
     * @param pos Where to add the atom, given in MolModel coordinates
     * @param bound_to The existing atom, if any, that the new atom should be
     * bound to.
     */
    virtual void addAtom(const RDGeom::Point3D& pos,
                         const RDKit::Atom* const bound_to = nullptr) = 0;

    /**
     * Add two new atoms and draw a bond between them
     * @param pos1 The location of the first atom in MolModel coordinates
     * @param pos2 The location of the second atom in MolModel coordinates
     */
    virtual void addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                                  const RDGeom::Point3D& pos2) = 0;

    /**
     * Add a bond between the two given atoms
     */
    virtual void addBond(const RDKit::Atom* const start_atom,
                         const RDKit::Atom* const end_atom) = 0;
};

} // namespace sketcher
} // namespace schrodinger
