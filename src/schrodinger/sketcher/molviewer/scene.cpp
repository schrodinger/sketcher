#include "schrodinger/sketcher/molviewer/scene.h"

#include <functional>
#include <unordered_set>

#include <GraphMol/CoordGen.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <QtGlobal>
#include <QApplication>
#include <QClipboard>
#include <QFont>
#include <QGraphicsSceneMouseEvent>
#include <QMimeData>
#include <QPainterPath>
#include <QString>
#include <QTransform>
#include <QUndoStack>
#include <QUrl>

#include "schrodinger/rdkit_extensions/convert.h"

#include "schrodinger/sketcher/file_import_export.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/qt_utils.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/dialog/file_save_image_dialog.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/stereochemistry.h"

using schrodinger::rdkit_extensions::Format;

#define SETTER_AND_GETTER(settings_member, update_method, type, getter, \
                          setter, variable_name)                        \
    type Scene::getter() const                                          \
    {                                                                   \
        return settings_member.variable_name;                           \
    }                                                                   \
    void Scene::setter(type value)                                      \
    {                                                                   \
        settings_member.variable_name = value;                          \
        update_method();                                                \
    }
#define ATOM_SETTING(type, getter, setter, variable_name)                 \
    SETTER_AND_GETTER(m_atom_item_settings, updateAtomAndBondItems, type, \
                      getter, setter, variable_name)
#define BOND_SETTING(type, getter, setter, variable_name)                  \
    SETTER_AND_GETTER(m_bond_item_settings, updateBondItems, type, getter, \
                      setter, variable_name)

namespace schrodinger
{
namespace sketcher
{

Scene::Scene(QObject* parent) : QGraphicsScene(parent)
{
    m_selection_highlighting_item = new SelectionHighlightingItem();
    m_predictive_highlighting_item = new PredictiveHighlightingItem();
    m_rect_select_item = new RectSelectionItem();
    m_ellipse_select_item = new EllipseSelectionItem();
    m_lasso_select_item = new LassoSelectionItem();
    m_undo_stack = new QUndoStack(this);
    m_mol_model = new MolModel(m_undo_stack, this);

    addItem(m_selection_highlighting_item);
    addItem(m_predictive_highlighting_item);
    // We intentionally don't add the m_*_select_item items to the scene.  Those
    // only get added during a drag select.

    connect(m_mol_model, &MolModel::moleculeChanged, this,
            &Scene::updateInteractiveItems);
    connect(m_mol_model, &MolModel::selectionChanged, this,
            &Scene::onMolModelSelectionChanged);
}

Scene::~Scene()
{
    // Avoid selection update calls when any selected item is destroyed
    disconnect(this, &Scene::selectionChanged, this,
               &Scene::updateSelectionHighlighting);
    // The QGraphicsScene destructor will destroy m_selection_highlighting_item
    // and m_predictive_highlighting_item for us, since those items are in the
    // scene
    // The m_*_select_item items shouldn't be in the scene unless this instance
    // is getting destroyed in the middle of a click-and-drag, but double check
    // just to be safe and make sure we avoid a double free.
    removeSelectItemsFromScene();
    delete m_rect_select_item;
    delete m_ellipse_select_item;
    delete m_lasso_select_item;
}

void Scene::removeSelectItemsFromScene()
{
    if (m_rect_select_item->scene() == this) {
        removeItem(m_rect_select_item);
    }
    if (m_ellipse_select_item->scene() == this) {
        removeItem(m_ellipse_select_item);
    }
    if (m_lasso_select_item->scene() == this) {
        removeItem(m_lasso_select_item);
    }
}

void Scene::setModel(SketcherModel* model)
{
    if (m_sketcher_model != nullptr) {
        throw std::runtime_error("The model has already been set");
    }
    m_sketcher_model = model;

    // Connect content-based signals
    connect(this, &Scene::changed, m_sketcher_model,
            &SketcherModel::interactiveItemsChanged);
    connect(m_sketcher_model, &SketcherModel::interactiveItemsRequested, this,
            [this]() { return getInteractiveItems(); });

    // Connect selection-based signals
    connect(this, &Scene::selectionChanged, m_sketcher_model,
            &SketcherModel::selectionChanged);
    connect(m_sketcher_model, &SketcherModel::selectionRequested, this,
            [this]() { return selectedItems(); });

    // Connect context menu-based signals
    connect(m_sketcher_model, &SketcherModel::contextMenuObjectsRequested, this,
            [this]() { return m_context_menu_objects; });
}

void Scene::loadMol(const RDKit::ROMol& mol)
{
    auto shared_mol = std::make_shared<RDKit::ROMol>(mol);
    loadMol(shared_mol);
}

void Scene::loadMol(std::shared_ptr<RDKit::ROMol> mol)
{
    m_mol_model->addMol(*mol);
}

void Scene::updateInteractiveItems()
{
    clearAllInteractiveItems();
    const auto* mol = m_mol_model->getMol();

    m_atom_to_atom_item.clear();
    m_bond_to_bond_item.clear();

    // create atom items
    for (std::size_t i = 0; i < mol->getNumAtoms(); ++i) {
        const auto* atom = mol->getAtomWithIdx(i);
        const auto pos = mol->getConformer().getAtomPos(i);
        auto* atom_item = new AtomItem(atom, m_fonts, m_atom_item_settings);
        atom_item->setPos(to_scene_xy(pos));
        m_atom_to_atom_item[atom] = atom_item;
        addItem(atom_item);
    }

    // create bond items
    for (auto bond : mol->bonds()) {
        const auto* from_atom_item = m_atom_to_atom_item[bond->getBeginAtom()];
        const auto* to_atom_item = m_atom_to_atom_item[bond->getEndAtom()];
        auto* bond_item = new BondItem(bond, *from_atom_item, *to_atom_item,
                                       m_fonts, m_bond_item_settings);
        m_bond_to_bond_item[bond] = bond_item;
        addItem(bond_item);
    }
}

std::shared_ptr<RDKit::ROMol> Scene::getRDKitMolecule() const
{
    // note that this will return a copy
    return std::make_shared<RDKit::ROMol>(*m_mol_model->getMol());
}

void Scene::importText(const std::string& text, Format format)
{
    auto show_import_failure = [&](const auto& exc) {
        auto text = QString("Import Failed: ") + exc.what();
        show_error_dialog("Import Error", text, window());
    };

    boost::shared_ptr<RDKit::RWMol> mol{nullptr};
    try {
        mol = rdkit_extensions::to_rdkit(text, format);
    } catch (const std::invalid_argument& exc) {
        show_import_failure(exc);
        return;
    } catch (const std::runtime_error& exc) {
        show_import_failure(exc);
        return;
    }

    // TODO: deal with chiral flag viz. SHARED-8774
    // TODO: honor existing coordinates if present
    RDKit::CoordGen::addCoords(*mol);

    assign_CIP_labels(*mol);

    loadMol(*mol);
}

std::string Scene::exportText(Format format)
{
    // TODO: handle reactions
    // TODO: handle toggling between selection/everything
    // TODO: handle export of selection as atom/bond properties
    return rdkit_extensions::to_string(*m_mol_model->getMol(), format);
}

// TODO: remove this method as part of SKETCH-1947
void Scene::clearInteractiveItems()
{
    m_mol_model->clear();
}

QList<QGraphicsItem*> Scene::getInteractiveItems() const
{
    // We expect all objects that inherit from AbstractGraphicsItem to
    // be interactive (atoms, bonds, etc.), and all objects that don't
    // to be purely graphical (selection highlighting paths, etc.)
    QList<QGraphicsItem*> interactive_items;
    for (auto item : items()) {
        if (dynamic_cast<AbstractGraphicsItem*>(item) != nullptr) {
            interactive_items.append(item);
        }
    }
    return interactive_items;
}

void Scene::clearAllInteractiveItems()
{
    // reset the state of the selection item
    removeSelectItemsFromScene();
    m_mouse_drag_action = MouseDragAction::NONE;
    m_lasso_select_item->clearPath();

    // remove all interactive items and reset the rdkit molecule; this will
    // preserve items include selection paths, highlighting items, and
    // potentially the watermark managed by the SketcherWidget
    for (auto item : getInteractiveItems()) {
        removeItem(item);
        delete item;
    }
}

// TODO: remove this method as part of SKETCH-1947
void Scene::selectAll()
{
    m_mol_model->selectAll();
}

// TODO: remove this method as part of SKETCH-1947
void Scene::invertSelection()
{
    m_mol_model->invertSelection();
}

// TODO: remove this method as part of SKETCH-1947
void Scene::clearSelectionPublic()
{
    m_mol_model->clearSelection();
}

void Scene::onImportTextRequested(const std::string& text, Format format)
{
    if (m_sketcher_model->getValueBool(
            ModelKey::NEW_STRUCTURES_REPLACE_CONTENT)) {
        clearInteractiveItems();
    }
    importText(text, format);
}

void Scene::showFileExportDialog()
{
    auto dialog = new FileExportDialog(m_sketcher_model, window());
    connect(dialog, &FileExportDialog::exportTextRequested, this,
            [this](Format format) {
                return QString::fromStdString(exportText(format));
            });
    dialog->show();
}

void Scene::showFileSaveImageDialog()
{
    auto dialog = new FileSaveImageDialog(m_sketcher_model, window());
    connect(dialog, &FileSaveImageDialog::exportImageRequested, this,
            [this](auto format, const auto& opts) {
                // TODO: this call uses the existing sketcherScene class;
                // update to render directly from this Scene instance
                return get_image_bytes(*m_mol_model->getMol(), format, opts);
            });
    dialog->show();
}

void Scene::onPasteRequested()
{
    auto data = QApplication::clipboard()->mimeData();
    if (data->hasText()) {
        importText(data->text().toStdString(), Format::AUTO_DETECT);
    }
}

qreal Scene::fontSize() const
{
    return m_fonts.size();
}

void Scene::setFontSize(qreal size)
{
    m_fonts.setSize(size);
    updateAtomAndBondItems();
}

ATOM_SETTING(CarbonLabels, carbonsLabeled, setCarbonsLabeled, m_carbon_labels)
ATOM_SETTING(bool, valenceErrorsShown, setValenceErrorsShown,
             m_valence_errors_shown)
BOND_SETTING(qreal, bondWidth, setBondWidth, m_bond_width)
BOND_SETTING(qreal, doubleBondSpacing, setDoubleBondSpacing,
             m_double_bond_spacing)

void Scene::updateAtomAndBondItems()
{
    // we need to update all of the atom items before we update any bond items,
    // since bond items pull information from their associated atom items
    for (auto item : items()) {
        if (auto atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            atom_item->updateCachedData();
        }
    }
    updateBondItems();
}

void Scene::updateBondItems()
{
    for (auto item : items()) {
        if (auto bond_item = qgraphicsitem_cast<BondItem*>(item)) {
            bond_item->updateCachedData();
        }
    }
}

void Scene::updateSelectionHighlighting()
{
    auto get_sel_path = [](AbstractGraphicsItem* item) {
        return item->selectionHighlightingPath();
    };
    QPainterPath path =
        buildHighlightingPathForItems(selectedItems(), get_sel_path);
    m_selection_highlighting_item->setHighlightingPath(path);
}

QPainterPath Scene::buildHighlightingPathForItems(
    QList<QGraphicsItem*> items,
    std::function<QPainterPath(AbstractGraphicsItem*)> path_getter) const
{
    QPainterPath path;
    for (auto item : items) {
        if (auto* molviewer_item = dynamic_cast<AbstractGraphicsItem*>(item)) {
            QPainterPath local_path = path_getter(molviewer_item);
            path |= molviewer_item->mapToScene(local_path);
        }
    }
    return path;
}

void Scene::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    m_mouse_down_scene_pos = event->scenePos();
    m_mouse_down_screen_pos = event->screenPos();
    if (event->button() == Qt::LeftButton &&
        m_sketcher_model->getDrawTool() == DrawTool::SELECT &&
        m_sketcher_model->getSelectionTool() == SelectionTool::LASSO) {
        m_lasso_select_item->clearPath();
        m_lasso_select_item->addPoint(m_mouse_down_scene_pos);
    }
}

void Scene::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    //  TODO: SKETCH-1890: translate on right click and rotate on middle click
    auto draw_tool = m_sketcher_model->getDrawTool();
    if (draw_tool == DrawTool::SELECT && (event->buttons() & Qt::LeftButton)) {
        auto select_tool = m_sketcher_model->getSelectionTool();
        // update the lasso path, even if the mouse hasn't yet moved far enough
        // to start a selection
        if (select_tool == SelectionTool::LASSO) {
            m_lasso_select_item->addPoint(event->scenePos());
        }

        // check to see whether the mouse has moved far enough to start a
        // selection
        if (m_mouse_drag_action == MouseDragAction::NONE) {
            int drag_dist = (m_mouse_down_screen_pos - event->screenPos())
                                .manhattanLength();
            if (drag_dist >= QApplication::startDragDistance()) {
                if (select_tool == SelectionTool::RECTANGLE) {
                    m_mouse_drag_action = MouseDragAction::RECTANGLE_SELECT;
                    addItem(m_rect_select_item);
                } else if (select_tool == SelectionTool::ELLIPSE) {
                    m_mouse_drag_action = MouseDragAction::ELLIPSE_SELECT;
                    addItem(m_ellipse_select_item);
                } else if (select_tool == SelectionTool::LASSO) {
                    m_mouse_drag_action = MouseDragAction::LASSO_SELECT;
                    addItem(m_lasso_select_item);
                }
            }
        }

        // update the marquee selection
        QRectF rect =
            rect_for_points(m_mouse_down_scene_pos, event->scenePos());
        if (m_mouse_drag_action == MouseDragAction::RECTANGLE_SELECT) {
            m_rect_select_item->setRect(rect);
            setPredictiveHighlightingForSelection(m_rect_select_item);
        } else if (m_mouse_drag_action == MouseDragAction::ELLIPSE_SELECT) {
            m_ellipse_select_item->setRect(rect);
            setPredictiveHighlightingForSelection(m_ellipse_select_item);
        } else if (m_mouse_drag_action == MouseDragAction::LASSO_SELECT) {
            setPredictiveHighlightingForSelection(m_lasso_select_item);
        }
    } else {
        // update predictive highlighting for mouse overs
        setPredictiveHighlightingForPoint(event->scenePos());
    }
}

void Scene::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    auto draw_tool = m_sketcher_model->getDrawTool();
    if (draw_tool == DrawTool::SELECT && event->button() == Qt::LeftButton) {
        SelectMode select_mode;
        auto modifiers = event->modifiers();
        if (modifiers & Qt::ControlModifier) {
            select_mode = SelectMode::TOGGLE;
        } else if (modifiers & Qt::ShiftModifier) {
            select_mode = SelectMode::SELECT;
        } else {
            select_mode = SelectMode::SELECT_ONLY;
        }

        QList<QGraphicsItem*> items_to_select;
        if (m_mouse_drag_action == MouseDragAction::NONE) {
            // this was a click, not a drag
            QGraphicsItem* item = itemAt(event->scenePos(), QTransform());
            if (item) {
                items_to_select.append(item);
            }
        } else {
            QAbstractGraphicsShapeItem* select_item;
            if (m_mouse_drag_action == MouseDragAction::RECTANGLE_SELECT) {
                select_item = m_rect_select_item;
            } else if (m_mouse_drag_action == MouseDragAction::ELLIPSE_SELECT) {
                select_item = m_ellipse_select_item;
            } else { // m_mouse_drag_action == MouseDragAction::LASSO_SELECT
                select_item = m_lasso_select_item;
            }

            // select the items within the marquee or lasso
            items_to_select = itemsWithinSelection(select_item);

            // reset the scene state
            removeItem(select_item);
            m_mouse_drag_action = MouseDragAction::NONE;
            m_predictive_highlighting_item->clearHighlightingPath();
            m_lasso_select_item->clearPath();
        }
        selectGraphicsItems(items_to_select, select_mode);
    }

    m_mouse_down_scene_pos = QPointF();
    m_mouse_down_screen_pos = QPointF();
}

void Scene::setPredictiveHighlightingForPoint(const QPointF& scene_pos)
{
    // TODO: highlight the context menu items if a context menu is visible
    // TODO: if the rotation tool is selected, highlight the entire selection
    //       when inside the rotation box
    QGraphicsItem* item = itemAt(scene_pos, QTransform());
    if (auto* highlightable_item = dynamic_cast<AbstractGraphicsItem*>(item)) {
        QPainterPath path = highlightable_item->predictiveHighlightingPath();
        path = highlightable_item->mapToScene(path);
        // note that setHighlightingPath is a no-op if path == the current
        // predictive highlighting path
        m_predictive_highlighting_item->setHighlightingPath(path);
    } else {
        m_predictive_highlighting_item->clearHighlightingPath();
    }
}

void Scene::setPredictiveHighlightingForSelection(
    const QAbstractGraphicsShapeItem* sel_item)
{
    QList<QGraphicsItem*> items_to_highlight = itemsWithinSelection(sel_item);
    auto get_pred_path = [](AbstractGraphicsItem* item) {
        return item->predictiveHighlightingPath();
    };
    QPainterPath path =
        buildHighlightingPathForItems(items_to_highlight, get_pred_path);
    m_predictive_highlighting_item->setHighlightingPath(path);
}

void Scene::selectGraphicsItems(const QList<QGraphicsItem*>& items,
                                const SelectMode select_mode)
{
    std::unordered_set<const RDKit::Atom*> atoms;
    std::unordered_set<const RDKit::Bond*> bonds;
    for (auto cur_item : items) {
        if (auto atom_item = qgraphicsitem_cast<AtomItem*>(cur_item)) {
            atoms.insert(atom_item->getAtom());
        } else if (auto bond_item = qgraphicsitem_cast<BondItem*>(cur_item)) {
            bonds.insert(bond_item->getBond());
        }
    }
    m_mol_model->select(atoms, bonds, select_mode);
}

void Scene::onMolModelSelectionChanged()
{
    // block the per-item signals for performance reasons
    Scene::SelectionChangeSignalBlocker signal_blocker(this);
    clearSelection();
    for (auto* atom : m_mol_model->getSelectedAtoms()) {
        m_atom_to_atom_item[atom]->setSelected(true);
    }
    for (auto* bond : m_mol_model->getSelectedBonds()) {
        m_bond_to_bond_item[bond]->setSelected(true);
    }
    updateSelectionHighlighting();
    // signal_blocker emits selectionChanged on destruction
}

QList<QGraphicsItem*>
Scene::itemsWithinSelection(const QAbstractGraphicsShapeItem* sel_item) const
{
    QList<QGraphicsItem*> items_within;
    for (auto item : items()) {
        if (sel_item->collidesWithItem(item)) {
            items_within.append(item);
        }
    }
    return items_within;
}

Scene::SelectionChangeSignalBlocker::SelectionChangeSignalBlocker(
    Scene* scene) :
    m_scene(scene),
    m_original_value(m_scene->signalsBlocked())
{
    m_scene->blockSignals(true);
}
Scene::SelectionChangeSignalBlocker::~SelectionChangeSignalBlocker()
{
    m_scene->blockSignals(m_original_value);
    m_scene->selectionChanged();
}

QWidget* Scene::window() const
{
    auto parent_widget = dynamic_cast<QWidget*>(parent());
    return parent_widget ? parent_widget->window() : nullptr;
}

void Scene::dropEvent(QGraphicsSceneDragDropEvent* event)
{
    auto mime_data = event->mimeData();
    for (auto url : mime_data->urls()) {
        if (url.isLocalFile()) { // File drop event
            auto file_path = url.toLocalFile();
            auto contents = get_file_text(file_path.toStdString());
            auto format = get_file_format(file_path);
            importText(contents, format);
        }
    }
    event->acceptProposedAction();
}

void Scene::dragEnterEvent(QGraphicsSceneDragDropEvent* event)
{
    event->acceptProposedAction();
}

void Scene::dragMoveEvent(QGraphicsSceneDragDropEvent* event)
{
    event->acceptProposedAction();
}

void Scene::dragLeaveEvent(QGraphicsSceneDragDropEvent* event)
{
    event->acceptProposedAction();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/scene.moc"
