/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <string>

#include <QtGlobal>
#include <QGraphicsScene>
#include <QPolygonF>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/molviewer/predictive_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/selection_items.h"

class QObject;
class QFont;

namespace RDKit
{
class ROMol;
}

namespace schrodinger
{

namespace rdkit_extensions
{
enum class Format;
}

namespace sketcher
{

class AbstractGraphicsItem;
class SketcherModel;
enum class DrawTool;
enum class ImageFormat;
enum class SelectionTool;
struct RenderOptions;

enum class MouseDragAction {
    NONE,
    RECTANGLE_SELECT,
    ELLIPSE_SELECT,
    LASSO_SELECT,
};

/**
 * A Qt graphics scene for displaying molecules.
 */
class SKETCHER_API Scene : public QGraphicsScene
{
    Q_OBJECT
  public:
    Scene(QObject* parent = nullptr);
    virtual ~Scene();

    /**
     * @param model A model instance to assign to this view
     * NOTE: the model must be set exactly once.
     */
    void setModel(SketcherModel* model);

    /**
     * Load an RDKit molecule into this scene
     *
     * @param mol The RDKit molecule to load.  Note that this class will store a
     * copy of this molecule and will not take ownership of the passed-in
     * instance.
     */
    void loadMol(const RDKit::ROMol& mol);
    /**
     * Load an RDKit molecule into this scene
     *
     * @param mol A shared pointer to the RDKit molecule to load.  After this
     * method is called, the molecule (including the conformation) must no
     * longer be modified by the calling code.
     */
    void loadMol(std::shared_ptr<RDKit::ROMol> mol);

    /**
     * @return RDKit molecule underlying what is currently rendered
     */
    std::shared_ptr<RDKit::ROMol> getRDKitMolecule() const;

    /**
     * Create a molecule from a text string and load that into the scene.
     * Atomic coordinates will be automatically generated using coordgen.
     *
     * @param text input data to load
     * @param format format to parse
     */
    void importText(const std::string& text, rdkit_extensions::Format format);

    /**
     * @param format format to convert to
     * @return serialized representation of the sketcher contents
     */
    std::string exportText(rdkit_extensions::Format format);

    /**
     *  Clear the contents of the scene
     */
    void clear();

    /**
     * Import the given text into the scene; optionally clear beforehand
     *
     * @param text input data to load
     * @param format format to parse
     */
    void onImportTextRequested(const std::string& text,
                               rdkit_extensions::Format format);

    /**
     * Present the user with a "Save Image" dialog.
     */
    void showFileSaveImageDialog();

    /**
     * Present the user with an "Export to File" dialog.
     */
    void showFileExportDialog();

    /**
     * Paste clipboard content into the scene
     */
    void onPasteRequested();

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

  protected:
    /**
     * If the m_rect_select_item or m_ellipse_select_item graphics items are in
     * the scene, remove them.  These items are normally removed in
     * mouseReleaseEvent, but this method can be used to ensure that they are
     * removed prior to a QGraphicsScene::clear() call or object destruction.
     */
    void removeSelectItemsFromScene();

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

    /**
     * Build a painter path for use with either predictive or selection
     * highlighting.
     *
     * @param items The graphics items that the path should include
     * @param path_getter A function to fetch the path for a given graphics item
     * @return The newly constructed path
     */
    QPainterPath buildHighlightingPathForItems(
        QList<QGraphicsItem*> items,
        std::function<QPainterPath(AbstractGraphicsItem*)> path_getter) const;

    // Override the QGraphicsScene mouse event methods
    void mousePressEvent(QGraphicsSceneMouseEvent* mouseEvent) override;
    void mouseMoveEvent(QGraphicsSceneMouseEvent* event) override;
    void mouseReleaseEvent(QGraphicsSceneMouseEvent* event) override;

    /**
     * Update predictive highlighting to highlight the item at the specified
     * point.  If the point is not over an item, predictive highlighting will be
     * cleared.
     */
    void setPredictiveHighlightingForPoint(const QPointF& scene_pos);

    /**
     * Update predictive highlighting to highlight all items within the user's
     * current marquee or lasso selection.
     *
     * @param sel_item The graphics item representing the current selection.
     */
    void setPredictiveHighlightingForSelection(
        const QAbstractGraphicsShapeItem* sel_item);

    /**
     * Find all graphics item within the user's current marquee or lasso
     * selection.
     *
     * @param sel_item The graphics item representing the current selection.
     */
    QList<QGraphicsItem*>
    itemsWithinSelection(const QAbstractGraphicsShapeItem* sel_item) const;

    std::shared_ptr<RDKit::ROMol> m_mol;
    Fonts m_fonts;
    AtomItemSettings m_atom_item_settings;
    BondItemSettings m_bond_item_settings;
    SketcherModel* m_sketcher_model = nullptr;
    SelectionHighlightingItem* m_selection_highlighting_item = nullptr;
    PredictiveHighlightingItem* m_predictive_highlighting_item = nullptr;
    RectSelectionItem* m_rect_select_item = nullptr;
    EllipseSelectionItem* m_ellipse_select_item = nullptr;
    LassoSelectionItem* m_lasso_select_item = nullptr;
    QPointF m_mouse_down_scene_pos;
    QPointF m_mouse_down_screen_pos;
    MouseDragAction m_mouse_drag_action = MouseDragAction::NONE;

    /**
     * Objects associated with the context menu instance that is currently open.
     * If no context menu instance is open (or if it is not associated with any
     * objects) then this should be empty.
     */
    QList<QGraphicsItem*> m_context_menu_objects;

  private:
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
