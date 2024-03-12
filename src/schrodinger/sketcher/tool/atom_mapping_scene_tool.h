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
};

} // namespace sketcher
} // namespace schrodinger
