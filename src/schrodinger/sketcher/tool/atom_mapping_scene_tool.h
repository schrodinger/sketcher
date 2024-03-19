#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"
#include <QGraphicsLineItem>
#include <QGraphicsPolygonItem>

namespace schrodinger
{
namespace sketcher
{

class Scene;
class MolModel;
class AtomItem;

class ArrowHeadItem : public QGraphicsPolygonItem
{
  public:
    ArrowHeadItem();
};

class ArrowLineItem : public QGraphicsLineItem
{
  public:
    ArrowLineItem();
};

enum class MappingAction : bool { ADD, REMOVE };
/**
 * A scene tool for adding and removing atom mapping numbers
 */
class SKETCHER_API AtomMappingSceneTool : public StandardSceneToolBase
{
  public:
    AtomMappingSceneTool(const MappingAction& action, Scene* scene,
                         MolModel* mol_model);

    // reimplemented AbstractSceneTool method
    void onLeftButtonPress(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    QPixmap createDefaultCursorPixmap() const override;

    /**
     * update the arrow items
     * @param line the segment to use for the arrow line. Arrow head will be at
     * line.p2()
     */
    void updateArrowItems(const QLineF& line);

    /**
     * @return the lowest number that is not being used as an atom mapping
     * number
     */
    int findLowestAvailableMappingNumber() const;

  protected:
    ArrowHeadItem m_arrow_head_item = ArrowHeadItem();
    ArrowLineItem m_arrow_line_item = ArrowLineItem();

    AtomItem* m_pressed_atom_item = nullptr;
    AtomItem* m_release_atom_item = nullptr;
    MappingAction m_action;

    /**
     * Determine whether the click-and-drag atom items represent a pair of atoms
     * that can be mapped. This will be true if and only if:
     *   - the hovered atom item is not nullptr
     *   - the two atom items are distinct (i.e. the user isn't hovering where
     *     the drag was started)
     *   - one of the atom items represents a reactant atom and one of the atom
     *     items represents a product atom.
     * @param pressed_atom_item The atom item where the click-and-drag was
     * started. This atom item must not be nullptr.
     * @param hovered_atom_item The atom item where the cursor currently is.
     * This atom item may be nullptr if the cursor isn't over an atom.
     */
    bool isValidMappingPair(const AtomItem* const pressed_atom_item,
                            const AtomItem* const hovered_atom_item);
};

} // namespace sketcher
} // namespace schrodinger
