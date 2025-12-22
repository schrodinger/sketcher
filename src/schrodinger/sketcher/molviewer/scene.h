#pragma once

#include <memory>
#include <string>
#include <unordered_map>
#include <variant>

#include <QGraphicsScene>
#include <QPolygonF>
#include <QtGlobal>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_items.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

class QObject;
class QFont;
class QPointF;

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

class AbstractGraphicsItem;
class AtomItem;
class BondItem;
class MolModel;
class NonMolecularItem;
class RotationItem;
class SGroupItem;
class SketcherModel;
enum class DrawTool;
enum class ImageFormat;
enum class ModelKey;
enum class SelectionTool;
struct RenderOptions;

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
     * @param selection_only If true, only consider selected items
     */
    QRectF getInteractiveItemsBoundingRect(
        const InteractiveItemFlagType types = InteractiveItemFlag::ALL,
        bool selection_only = false) const;

    /**
     * @return the bounding rectangle of all items in the scene. This includes
     * all interactive items plus any annotation items such as the simplified
     * stereo label
     */
    QRectF getSceneItemsBoundingRect() const;

    /**
     * get atoms, bonds, secondary connections, sgroups, and non-molecular
     * objects in the scene
     * @param subset The subset of objects to return (all, selected,
     * selected_or_hovered).  Defaults to returning all objects.
     * @param pos The position of the mouse cursor, this is needed to figure
     * out which items are hovered
     * @return a tuple of sets of atoms, bonds, secondary connections, sgroups,
     * and non-molecular objects
     */
    std::tuple<std::unordered_set<const RDKit::Atom*>,
               std::unordered_set<const RDKit::Bond*>,
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

    //*
    /* @return a list of all graphics items in the scene that collied with the
     * given item, excluding bond items whose center doesn't fall within the
     * bounds of the item
     */
    QList<QGraphicsItem*>
    getCollidingItemsUsingBondMidpoints(QGraphicsItem* item) const;

    /**
     * Make sure that the returned list includes graphics items for both an
     * attachment point dummy atom *and* the associated bond if either the atom
     * *or* the bond are included in the input list.
     */
    QList<QGraphicsItem*>
    ensureCompleteAttachmentPoints(const QList<QGraphicsItem*>& items) const;

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

    /**
     * @return true if the user is in the middle of a drag-and-drop rotation or
     * translation.  False otherwise.
     */
    bool isDuringAtomDrag();

    /**
     * Update the graphical items colors in response to a change in the
     * background color (e.g. switching to/from dark mode).
     */

    void onBackgroundColorChanged();

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
        const std::unordered_set<const RDKit::Bond*>& secondary_connections,
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
        const std::unordered_set<const NonMolecularObject*>&
            non_molecular_objects);

    void viewportTranslationRequested(const QPointF& start_screen_position,
                                      const QPointF& end_screen_position);

    /**
     * Notify the view of the appropriate cursor hint for the current
     * left-button scene tool
     */
    void newCursorHintRequested(const QPixmap& pixmap);

    /**
     * Emitted when a drag-and-drop operation has completed *if and only if* the
     * drag-and-drop did not affect atom connectivity (i.e. if no atoms were
     * merged as a result of the rotation or translation).  If connectivity was
     * affected, then MolModel::modelChanged will notify listeners of the
     * change.
     */
    void representationChangingAtomDragFinished();

    void atomHovered(const RDKit::Atom*);
    void bondHovered(const RDKit::Bond*);

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
     * Update all atoms and bond graphics items and/or arrow and plus graphics
     * items and regenerate them from the MolModel data
     * @param what_changed The type of change that occurred in the MolModel
     */
    void updateItems(const WhatChangedType what_changed);

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
     * in any way.). Also show or hide simplified stereochemistry annotations
     * based on the current settings.
     */
    void onDisplaySettingsChanged();

    /**
     * calculate a preview of monomer label sizes for all monomeric atoms and
     * store them as RDKit properties. This is needed in monomer coordgen to
     * allocate enough space for the labels and avoid overlaps
     */
    void updateMonomerLabelSizeOnModel();

    /**
     * Update the path drawn to show selection highlighting.
     */
    void updateSelectionHighlighting();

    /**
     * Update the current tool when the selection changes
     */
    void updateSceneToolAfterSelection();

    /**
     * Update the path drawn to show item colored highlights
     */
    void updateHaloHighlighting();

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

    void onAtomDragStarted();
    void onAtomDragFinished(const bool were_atoms_merged);

    void updateHovered(const QPointF& pos);
    void setHovered(
        std::variant<std::monostate, const RDKit::Atom*, const RDKit::Bond*>
            new_hovered);
    void clearHovered();

    /**
     * Set the scene tool (i.e. the mouse cursor mode) to the given value
     */
    void setSceneTool(std::shared_ptr<AbstractSceneTool> new_scene_tool);

    Fonts m_fonts;
    MolModel* m_mol_model = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
    SelectionHighlightingItem* m_selection_highlighting_item = nullptr;
    // a QGraphicsItem that groups all the user halo highlights. Each item can
    // highlight atoms and bonds with a single color, so we might need multiple
    // child items
    QGraphicsItemGroup* m_halo_highlighting_item = nullptr;
    QGraphicsTextItem* m_simplified_stereo_label = nullptr;
    QPointF m_mouse_down_screen_pos;

    /// A set of all "interactive" graphics items that are currently in the
    /// scene.  See the getInteractiveItems docstring for an explanation of
    /// "interactive" versus "non-interactive" items.
    std::unordered_set<QGraphicsItem*> m_interactive_items;

    std::unordered_map<const RDKit::Atom*, QGraphicsItem*> m_atom_to_atom_item;
    std::unordered_map<const RDKit::Bond*, QGraphicsItem*> m_bond_to_bond_item;
    /**
     * In monomeric models, some bonds represent more than one connection
     * between two monomers (e.g. neighboring cysteines additionally joined by a
     * disulfide bond). We keep track of graphics items that represent secondary
     * connections separately.
     */
    std::unordered_map<const RDKit::Bond*, QGraphicsItem*>
        m_bond_to_secondary_connection_item;
    std::unordered_map<const RDKit::SubstanceGroup*, SGroupItem*>
        m_s_group_to_s_group_item;
    std::unordered_map<const NonMolecularObject*, NonMolecularItem*>
        m_non_molecular_to_non_molecular_item;
    std::shared_ptr<AbstractSceneTool> m_scene_tool;
    /**
     * Some mouse events require m_scene_tool to change (e.g. clicking with the
     * add reaction arrow tool equips the add plus tool). In order for the tool
     * to be still present when the mouse is released, we store the tool that
     * was active when the mouse was pressed and clear it on mouse release.
     */
    std::shared_ptr<AbstractSceneTool> m_scene_tool_from_mouse_press;

    /**
     * The atom or bond (or nothing) that the cursor is currently over
     */
    std::variant<std::monostate, const RDKit::Atom*, const RDKit::Bond*>
        m_hovered = std::monostate();

    /**
     * Objects associated with the context menu instance that is currently open.
     * If no context menu instance is open (or if it is not associated with any
     * objects) then this should be empty.
     */
    std::vector<const RDKit::Atom*> m_context_menu_atoms;
    bool m_currently_dragging_atom;

    /**
     * Check if we are currently displaying a simplified stereo annotation
     */
    bool isSimplifiedStereoAnnotationVisible() const;

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

// This function is defined here rather than in image_generation.h along with
// the other overloads to avoid introducing a dependency on private classes
// (Scene) in the public header file. The implementation is in
// image_generation.cpp

/**
 * @param Scene to render
 * @param format format of the image
 * @param opts given image generation configuration
 * @return byte array of data generated from the 2D sketcher
 */
SKETCHER_API QByteArray get_image_bytes(Scene& scene, ImageFormat format,
                                        const RenderOptions& opts);

} // namespace sketcher
} // namespace schrodinger
