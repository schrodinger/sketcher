#pragma once

#include <tuple>
#include <utility>

#include <QGraphicsPathItem>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/abstract_draw_atom_bond_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A blue squiggle used to show where an attachment would be added if the user
 * clicked or released their mouse drag.
 */
class HintSquiggleItem : public QGraphicsPathItem
{
  public:
    HintSquiggleItem(QGraphicsItem* parent = nullptr);
};

/**
 * A scene tool for drawing attachment points
 */
class SKETCHER_API DrawAttachmentPointSceneTool : public AbstractDrawSceneTool
{
  public:
    DrawAttachmentPointSceneTool(Scene* scene, MolModel* mol_model);

  protected:
    HintSquiggleItem m_hint_squiggle_item = HintSquiggleItem();

    // reimplemented AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;

    QPixmap createDefaultCursorPixmap() const override;

    // reimplemented AbstractDrawSceneTool methods
    void onBondClicked(const RDKit::Bond* const bond) override;

    bool
    shouldDrawBondForClickOnAtom(const RDKit::Atom* const atom) const override;

    void addAtom(const RDGeom::Point3D& pos,
                 const RDKit::Atom* const bound_to) override;

    std::pair<QPointF, const RDKit::Atom*>
    getBondEndInMousedDirection(const QPointF& start_pos,
                                const RDKit::Atom* const start_atom,
                                const QPointF& mouse_pos) const override;

    std::tuple<RDGeom::Point3D, QPointF, const RDKit::Atom*>
    getDefaultBondPosition(const RDKit::Atom* const atom) const override;

    std::tuple<bool, QPointF, const RDKit::Atom*>
    getDragStartInfo() const override;

    void setHintBondVisible(bool visible) override;

    void updateHintBondPath(const QPointF& start, const QPointF& end) override;

    // AbstractDrawSceneTool methods that are not relevant to this tool
    void onEmptySpaceClicked(const RDGeom::Point3D& pos) override{};

    void addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                          const RDGeom::Point3D& pos2) override{};

    void addBond(const RDKit::Atom* const start_atom,
                 const RDKit::Atom* const end_atom) override{};
};

} // namespace sketcher
} // namespace schrodinger
