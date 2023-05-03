#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include <QGraphicsScene>
#include <QPolygonF>
#include <QtGlobal>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/mol_model.h"
#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/scene_tools/abstract_scene_tool.h"
#include "schrodinger/sketcher/molviewer/scene_tools/select_scene_tool.h"
#include "schrodinger/sketcher/molviewer/selection_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_items.h"

class QObject;
class QFont;
class QUndoStack;

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

class AbstractGraphicsItem;
class AtomItem;
class BondItem;
class MolModel;
class SketcherModel;
enum class DrawTool;
enum class ModelKey;
enum class SelectionTool;

/**
 * A scene tool (i.e. a mouse cursor mode) that does nothing and draws nothing.
 * Setting a NullSceneTool on the Scene will safely unload the existing scene
 * tool.
 */
class NullSceneTool : public AbstractSceneTool
{
  public:
    NullSceneTool();
};

/**
 * A Qt graphics scene for displaying molecules.
 */
class SKETCHER_API Scene : public QGraphicsScene
{
    Q_OBJECT
  public:
    Scene(MolModel* mol_model, SketcherModel* sketcher_model,
          QObject* parent = nullptr);
    virtual ~Scene();

    /**
     * @return interactive items in the scene; these are items that inherit
     * from AbstractGraphicsItem (atoms, bonds, etc.) as opposed to objects
     * that are purely graphical (selection highlighting paths, etc.)
     */
    QList<QGraphicsItem*> getInteractiveItems() const;

    /**
     * Update the MolModel selection for the atoms and bonds that correspond to
     * the specified graphics items.  This will trigger a call to
     * onMolModelSelectionChanged, which is responsible for actually selecting
     * the graphics items.  This method should always be used to select graphics
     * items to ensure that selection is kept in sync between MolModel and the
     * Scene.
     * @param items The graphics items to update the selection of
     * @param select_mode Whether to select, deselect, toggle selection, or
     * select-only (i.e. clear the selection and then select)
     */
    void selectGraphicsItems(const QList<QGraphicsItem*>& items,
                             const SelectMode select_mode);

    // Getters and setters for changing settings
    qreal fontSize() const;
    void setFontSize(qreal size);

    bool allAtomsShown() const;
    void setAllAtomsShown(bool value);

    CarbonLabels carbonsLabeled() const;
    void setCarbonsLabeled(CarbonLabels value);

    bool valenceErrorsShown() const;
    void setValenceErrorsShown(bool value);

    qreal bondWidth() const;
    void setBondWidth(qreal value);

    qreal doubleBondSpacing() const;
    void setDoubleBondSpacing(qreal value);

  signals:
    void importTextRequested(const std::string& text,
                             rdkit_extensions::Format format);

  protected:
    using QGraphicsScene::clear;
    using QGraphicsScene::clearSelection;

    /**
     * Delete all interactive graphics items in the scene; these are items that
     * inherit from AbstractGraphicsItem (atoms, bonds, etc.) as opposed to
     * objects that are purely graphical (selection highlighting paths, etc.).
     */
    void clearAllInteractiveItems();

    /**
     * Clear all interactive graphics items (e.g. atoms, bonds, etc.) and
     * regenerate them from the current MolModel molecule.
     */
    void updateInteractiveItems();

    /**
     * Call updateCachedData() on all AtomItems and BondItems in the scene.
     * (BondItems always need updating after their bound AtomItems are modified
     * in any way.)
     */
    void updateAtomAndBondItems();

    /**
     * Call updateCachedData() on all BondItems (but not AtomItems) in the
     * scene.
     */
    void updateBondItems();

    /**
     * Update the path drawn to show selection highlighting.
     */
    void updateSelectionHighlighting();

    // Override the QGraphicsScene mouse event methods
    void mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent) override;
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event) override;

    /**
     * Update the graphics items selection in response to a change in the
     * MolModel selection.
     */
    void onMolModelSelectionChanged();

    /**
     * Update the scene in response to changes in the SketcherModel
     * @param keys The keys that were changed
     */
    void onModelValuesChanged(const std::unordered_set<ModelKey>& keys);

    /**
     * Update the scene tool (i.e. the mouse cursor mode) based on the current
     * SketcherModel settings
     */
    void updateSceneTool();

    /**
     * Set the scene tool (i.e. the mouse cursor mode) to the given value
     */
    void setSceneTool(std::shared_ptr<AbstractSceneTool> new_scene_tool);

    Fonts m_fonts;
    AtomItemSettings m_atom_item_settings;
    BondItemSettings m_bond_item_settings;
    MolModel* m_mol_model = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
    SelectionHighlightingItem* m_selection_highlighting_item = nullptr;
    QPointF m_mouse_down_screen_pos;
    std::unordered_map<const RDKit::Atom*, AtomItem*> m_atom_to_atom_item;
    std::unordered_map<const RDKit::Bond*, BondItem*> m_bond_to_bond_item;
    std::shared_ptr<AbstractSceneTool> m_scene_tool;
    bool m_drag_started = false;

    /**
     * Objects associated with the context menu instance that is currently open.
     * If no context menu instance is open (or if it is not associated with any
     * objects) then this should be empty.
     */
    QList<QGraphicsItem*> m_context_menu_objects;

  private:
    /**
     * Signal blocker for selection based signals
     */
    class SelectionChangeSignalBlocker
    {
      public:
        SelectionChangeSignalBlocker(Scene* scene);
        ~SelectionChangeSignalBlocker();

      private:
        Scene* m_scene;
        bool m_original_value;
    };

    /**
     * @return window for the parent widget
     */
    QWidget* window() const;

    /**
     * Overrides the drag/drop events to import text files of supported formats
     */
    void dragEnterEvent(QGraphicsSceneDragDropEvent* event) override;
    void dragMoveEvent(QGraphicsSceneDragDropEvent* event) override;
    void dragLeaveEvent(QGraphicsSceneDragDropEvent* event) override;
    void dropEvent(QGraphicsSceneDragDropEvent* event) override;
};

} // namespace sketcher
} // namespace schrodinger
