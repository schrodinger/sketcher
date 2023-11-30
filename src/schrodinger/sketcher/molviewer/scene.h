#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include <QGraphicsScene>
#include <QPolygonF>
#include <QtGlobal>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_items.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

class QObject;
class QFont;

namespace RDKit
{
class Atom;
class Bond;
class ROMol;
class SubstanceGroup;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class RotationItem;
class AbstractGraphicsItem;
class AtomItem;
class BondItem;
class MolModel;
class NonMolecularItem;
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
          QWidget* parent = nullptr);
    virtual ~Scene();

    /**
     * Get a list of interactive items in the scene; these are items that
     * inherit from AbstractGraphicsItem (atoms, bonds, etc.) as opposed to
     * objects that are purely graphical (selection highlighting paths, etc.)
     *
     * @param types The subset of interactive items to return.  Defaults to
     * returning all interactive items.
     * @return all requested items
     */
    QList<QGraphicsItem*> getInteractiveItems(
        const InteractiveItemFlagType types = InteractiveItemFlag::ALL) const;

    /**
     * @return a set of all selected interactive items in the scene
     */
    std::unordered_set<QGraphicsItem*> getSelectedInteractiveItems() const;

    /**
     * Get the bounding rectangle of all interactive items in the scene of a
     * specified type
     * @param types The subset of interactive items to return.  Defaults to
     * returning all interactive items.
     */
    QRectF getInteractiveItemsBoundingRect(
        const InteractiveItemFlagType types = InteractiveItemFlag::ALL) const;

    /**
     * get atoms, bonds, sgroups and non-molecular objects in the scene
     * @param subset The subset of objects to return (all, selected, selected or
    hovered).  Defaults to returning all
     * objects.
     * @param pos The position of the mouse cursor, this is needed to figure
    out which items are hovered
     * @return a tuple of sets of atoms, bonds, sgroups and non-molecular
     * objects
     */
    std::tuple<std::unordered_set<const RDKit::Atom*>,
               std::unordered_set<const RDKit::Bond*>,
               std::unordered_set<const RDKit::SubstanceGroup*>,
               std::unordered_set<const NonMolecularObject*>>
    getModelObjects(SceneSubset subset = SceneSubset::ALL,
                    QPointF* pos = nullptr) const;
    /**
     * Get the topmost interactive item at a given position.  See the
     * getInteractiveItems docstring for an explanation of "interactive"
     * items.
     *
     * @param pos The Scene coordinates to use
     * @param types The subset of interactive item types to return.  Any
     * item types that are not specified by this flag will be ignored.
     * @return the top interactive graphics item of the requested type(s) at
     * the given position, or nullptr if none is found
     */
    AbstractGraphicsItem*
    getTopInteractiveItemAt(const QPointF& pos,
                            const InteractiveItemFlagType types) const;

    /**
     * Make sure that the returned list includes graphics items for both an
     * attachment point dummy atom *and* the associated bond if either the atom
     * *or* the bond are included in the input list.
     */
    QList<QGraphicsItem*>
    ensureCompleteAttachmentPoints(const QList<QGraphicsItem*>& items) const;

    // These setters are public because the corresponding setting doesn't yet
    // exist in SketcherModel, so they're called directly from SketcherWidget.
    // These setters should be made private as part of SKETCH-2066.
    void setFontSize(const qreal size);
    void setCarbonsLabeled(const CarbonLabels value);

    /**
     * display the appropriate context menu at the given position
     * @param event The mouse event that triggered the context menu
     */
    void showContextMenu(QGraphicsSceneMouseEvent* event);

    QRectF getSelectionRect() const;

    /**
     * Notify the active left-click scene tool that the mouse cursor has left
     * the view.
     */
    void onMouseLeave();

    /**
     * Emit a signal notifying the view of the current left-button scene tool's
     * cursor
     */
    void requestCursorHintUpdate();

    AtomItem* getAtomItemForAtom(const RDKit::Atom* atom) const;

  signals:
    /**
     * Request that the widget import the given text in the given format
     */
    void importTextRequested(const std::string& text,
                             rdkit_extensions::Format format);

    /**
     * Request the widget show the appropriate context menu at the current
     * position given the relevant atoms/bonds/sgroups
     */
    void showContextMenuRequested(
        QGraphicsSceneMouseEvent* event,
        const std::unordered_set<const RDKit::Atom*>& atoms,
        const std::unordered_set<const RDKit::Bond*>& bonds,
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups);

    void viewportTranslationRequested(const QPointF& start_screen_position,
                                      const QPointF& end_screen_position);

    /**
     * Notify the view of the appropriate cursor hint for the current
     * left-button scene tool
     */
    void newCursorHintRequested(const QPixmap& pixmap);

  protected:
    using QGraphicsScene::clear;
    using QGraphicsScene::clearSelection;

    /**
     * Deletes interactive graphics items in the scene; these are items that
     * inherit from AbstractGraphicsItem (atoms, bonds, pluses, etc.) as opposed
     * to objects that are purely graphical (selection highlighting paths,
     * etc.).
     * @param types The subset of interactive items to delete.  Defaults to
     * deleting all interactive items.  Note that this method does *not*
     * distinguish between attachment points and non-attachment points, so if
     * this flag has the ATOM_NOT_AP (or BOND_NOT_AP) bit set, it must also have
     * the ATTACHMENT_POINT (or ATTACHMENT_POINT_BOND) bit set.  (In other
     * words, don't try to delete *only* attachment points using this method.)
     */
    void clearInteractiveItems(
        const InteractiveItemFlagType types = InteractiveItemFlag::ALL);

    /**
     * Clear all atom and bond graphics items and regenerate them from the
     * MolModel molecule.
     */
    void updateMolecularItems();

    /**
     * Clear all arrow and plus graphics items and regenerate them from the
     * MolModel.
     */
    void updateNonMolecularItems();

    /**
     * update the positions of all interactive items in the scene
     */
    void moveInteractiveItems();

    /**
     * functions to select atoms, bonds and non molecular items in the scene
     * according to the selection info stored in m_mol_model.
     */
    void updateItemSelection();

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
    void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event) override;

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

    void setColorHeteroatoms(const bool color_heteroatoms);

    /**
     * Return the scene tool that should be used for a mouse event, based on the
     * button (left, middle or right) that was pressed
     * @param event The mouse event
     */
    std::shared_ptr<AbstractSceneTool>
    getSceneTool(QGraphicsSceneMouseEvent* const event);

    /**
     * Update the scene tool (i.e. the mouse cursor mode) based on the
     * current SketcherModel settings
     */
    void updateSceneTool();

    /**
     * Instantiate and return a new scene tool (i.e. the mouse cursor mode)
     * based on the current SketcherModel settings.
     */
    std::shared_ptr<AbstractSceneTool> getNewSceneTool();

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

    /// A set of all "interactive" graphics items that are currently in the
    /// scene.  See the getInteractiveItems docstring for an explanation of
    /// "interactive" versus "non-interactive" items.
    std::unordered_set<QGraphicsItem*> m_interactive_items;

    std::unordered_map<const RDKit::Atom*, AtomItem*> m_atom_to_atom_item;
    std::unordered_map<const RDKit::Bond*, BondItem*> m_bond_to_bond_item;
    std::unordered_map<const NonMolecularObject*, NonMolecularItem*>
        m_non_molecular_to_non_molecular_item;
    std::shared_ptr<AbstractSceneTool> m_left_button_scene_tool;
    std::shared_ptr<AbstractSceneTool> m_middle_button_scene_tool;
    std::shared_ptr<AbstractSceneTool> m_right_button_scene_tool;

    bool m_drag_started = false;
    bool m_is_during_double_click = false;

    /**
     * Objects associated with the context menu instance that is currently open.
     * If no context menu instance is open (or if it is not associated with any
     * objects) then this should be empty.
     */
    std::vector<const RDKit::Atom*> m_context_menu_atoms;

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
