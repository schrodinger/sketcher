#include "schrodinger/sketcher/sketcher_widget.h"

#include <QClipboard>
#include <QCursor>
#include <QGraphicsPixmapItem>
#include <QMimeData>
#include <QScreen>
#include <QWidget>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <boost/algorithm/string.hpp>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/val.h>
#endif

#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/dialog/file_save_image_dialog.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/menu/atom_context_menu.h"
#include "schrodinger/sketcher/menu/background_context_menu.h"
#include "schrodinger/sketcher/menu/bond_context_menu.h"
#include "schrodinger/sketcher/menu/selection_context_menu.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/view.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/molops.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/rdkit/molops.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

SketcherWidget::SketcherWidget(QWidget* parent) :
    QWidget(parent),
    m_undo_stack(new QUndoStack(this)),
    m_mol_model(new MolModel(m_undo_stack)),
    m_sketcher_model(new SketcherModel(this)),
    m_scene(new Scene(m_mol_model, m_sketcher_model, this))
{
    // The tools in ~Scene will access the underlying mol, so we need to
    // make sure the mol model still exists when the scene is destroyed.
    // This is controlled by the order in which parentship relationships
    // are defined.
    m_mol_model->setParent(this);

    m_ui.reset(new Ui::SketcherWidgetForm());
    m_ui->setupUi(this);

    m_ui->top_bar_wdg->setModel(m_sketcher_model);
    m_ui->side_bar_wdg->setModel(m_sketcher_model);

    m_ui->view->setScene(m_scene);
    m_ui->view->setMolModel(m_mol_model);
    connect(m_scene, &Scene::importTextRequested, this,
            &SketcherWidget::importText);
    connect(m_scene, &Scene::showContextMenuRequested, this,
            &SketcherWidget::showContextMenu);

    // Connect the scene to the view
    connect(m_scene, &Scene::viewportTranslationRequested, m_ui->view,
            &View::translateViewportFromScreenCoords);
    connect(m_scene, &Scene::newCursorHintRequested, m_ui->view,
            &View::onNewCursorHintRequested);

    // Connect the SketcherModel to the undo stack
    connect(m_sketcher_model, &SketcherModel::undoStackCanUndoRequested,
            m_undo_stack, &QUndoStack::canUndo);
    connect(m_sketcher_model, &SketcherModel::undoStackCanRedoRequested,
            m_undo_stack, &QUndoStack::canRedo);
    connect(m_undo_stack, &QUndoStack::canUndoChanged, m_sketcher_model,
            &SketcherModel::undoStackDataChanged);
    connect(m_undo_stack, &QUndoStack::canRedoChanged, m_sketcher_model,
            &SketcherModel::undoStackDataChanged);

    connect(m_mol_model, &MolModel::reactionArrowAdded, m_sketcher_model,
            &SketcherModel::onReactionArrowAdded);
    connect(m_sketcher_model, &SketcherModel::reactionCountRequested,
            m_mol_model, &MolModel::hasReactionArrow);
    connect(m_sketcher_model, &SketcherModel::valuePinged, this,
            &SketcherWidget::onModelValuePinged);

    connectTopBarSlots();
    connectSideBarSlots();

    // Create and connect the context menus
    m_atom_context_menu = new AtomContextMenu(m_sketcher_model, this);
    m_bond_context_menu = new BondContextMenu(this);
    m_selection_context_menu = new SelectionContextMenu(m_sketcher_model, this);
    m_background_context_menu =
        new BackgroundContextMenu(m_sketcher_model, this);
    connectContextMenu(*m_atom_context_menu);
    connectContextMenu(*m_bond_context_menu);
    connectContextMenu(*m_selection_context_menu);
    connectContextMenu(*m_background_context_menu);

    setStyleSheet(schrodinger::sketcher::SKETCHER_WIDGET_STYLE);

    // Set up the watermark
    m_watermark_item = new QGraphicsPixmapItem();
    m_watermark_item->setPixmap(QPixmap(":icons/2D-Sketcher-watermark.svg"));
    m_watermark_item->setFlag(QGraphicsItem::ItemIgnoresTransformations, true);
    m_scene->addItem(m_watermark_item);
    connect(m_scene, &Scene::changed, this, &SketcherWidget::updateWatermark);
    connect(m_ui->view, &View::resized, this, &SketcherWidget::updateWatermark);

    // use the custom cursor everywhere
    setCursor(
        QCursor(get_arrow_cursor_pixmap(), CURSOR_HOTSPOT_X, CURSOR_HOTSPOT_Y));
    // force the scene to update the view's cursor now that all of the signals
    // are connected
    m_scene->requestCursorHintUpdate();
}

SketcherWidget::~SketcherWidget() = default;

void SketcherWidget::connectTopBarSlots()
{
    connect(m_ui->top_bar_wdg, &SketcherTopBar::undoRequested, m_undo_stack,
            &QUndoStack::undo);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::redoRequested, m_undo_stack,
            &QUndoStack::redo);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::cleanupRequested, [this]() {
        m_mol_model->regenerateCoordinates();
        m_ui->view->fitToScreen();
    });
    connect(m_ui->top_bar_wdg, &SketcherTopBar::fitToScreenRequested,
            m_ui->view, &View::fitToScreen);

    // Connect "More Actions" menu
    connect(m_ui->top_bar_wdg, &SketcherTopBar::flipHorizontalRequested,
            m_mol_model, &MolModel::flipAllHorizontal);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::flipVerticalRequested,
            m_mol_model, &MolModel::flipAllVertical);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::addExplicitHydrogensRequested,
            m_mol_model,
            [this]() { m_mol_model->updateExplicitHs(ExplicitHActions::ADD); });
    connect(
        m_ui->top_bar_wdg, &SketcherTopBar::removeExplicitHydrogensRequested,
        m_mol_model,
        [this]() { m_mol_model->updateExplicitHs(ExplicitHActions::REMOVE); });
    connect(m_ui->top_bar_wdg, &SketcherTopBar::selectAllRequested, m_mol_model,
            &MolModel::selectAll);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::clearSelectionRequested,
            m_mol_model, &MolModel::clearSelection);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::invertSelectionRequested,
            m_mol_model, &MolModel::invertSelection);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::cutRequested, this,
            &SketcherWidget::cut);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::copyRequested, this,
            &SketcherWidget::copy);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::pasteRequested, this,
            &SketcherWidget::paste);

    // Clear/Import/Export
    connect(m_ui->top_bar_wdg, &SketcherTopBar::clearSketcherRequested,
            m_mol_model, &MolModel::clear);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::importTextRequested, this,
            &SketcherWidget::importText);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::saveImageRequested, this,
            &SketcherWidget::showFileSaveImageDialog);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::exportToFileRequested, this,
            &SketcherWidget::showFileExportDialog);
}

void SketcherWidget::connectSideBarSlots()
{
    // Connect "Select Options" widget
    connect(m_ui->side_bar_wdg, &SketcherSideBar::selectAllRequested,
            m_mol_model, &MolModel::selectAll);
    connect(m_ui->side_bar_wdg, &SketcherSideBar::clearSelectionRequested,
            m_mol_model, &MolModel::clearSelection);
    connect(m_ui->side_bar_wdg, &SketcherSideBar::invertSelectionRequested,
            m_mol_model, &MolModel::invertSelection);
}

void SketcherWidget::connectContextMenu(const ModifyAtomsMenu& menu)
{
    // TODO
}

void SketcherWidget::connectContextMenu(const ModifyBondsMenu& menu)
{
    // TODO
}

void SketcherWidget::connectContextMenu(const SelectionContextMenu& menu)
{
    // TODO

    connectContextMenu(*menu.m_modify_atoms_menu);
    connectContextMenu(*menu.m_modify_bonds_menu);
}

void SketcherWidget::connectContextMenu(const BackgroundContextMenu& menu)
{
    connect(&menu, &BackgroundContextMenu::saveImageRequested, this,
            &SketcherWidget::showFileSaveImageDialog);
    connect(&menu, &BackgroundContextMenu::exportToFileRequested, this,
            &SketcherWidget::showFileExportDialog);
    connect(&menu, &BackgroundContextMenu::undoRequested, m_undo_stack,
            &QUndoStack::undo);
    connect(&menu, &BackgroundContextMenu::redoRequested, m_undo_stack,
            &QUndoStack::redo);
    connect(&menu, &BackgroundContextMenu::flipHorizontalRequested, m_mol_model,
            &MolModel::flipAllHorizontal);
    connect(&menu, &BackgroundContextMenu::flipVerticalRequested, m_mol_model,
            &MolModel::flipAllVertical);
    connect(&menu, &BackgroundContextMenu::selectAllRequested, m_mol_model,
            &MolModel::selectAll);
    connect(&menu, &BackgroundContextMenu::copyRequested, this,
            &SketcherWidget::copy);
    connect(&menu, &BackgroundContextMenu::pasteRequested, this,
            &SketcherWidget::paste);
    connect(&menu, &BackgroundContextMenu::clearRequested, m_mol_model,
            &MolModel::clear);
}

/**
 * @internal
 * @param mol_model the model to extract the molecule from
 * @param subset whether to extract everything or just the current selection
 */
static boost::shared_ptr<RDKit::ROMol> extract_mol(MolModel* mol_model,
                                                   SceneSubset subset)
{
    // TODO: strip TAG_PROPERTY from atoms/bonds before returning
    if (subset == SceneSubset::SELECTION) {
        return mol_model->getSelectedMol();
    }
    return boost::make_shared<RDKit::ROMol>(*mol_model->getMol());
}

/**
 * @internal
 * @param mol_model the model to extract the molecule from
 * @param format the format to serialize the molecule to
 * @param subset whether to extract everything or just the current selection
 */
static std::string extract_string(MolModel* mol_model, Format format,
                                  SceneSubset subset)
{
    if (format == Format::MDL_MOLV2000) {
        throw std::invalid_argument(
            "Sketcher does not support exporting MDL_MOLV2000 format");
    }

    if (mol_model->hasReactionArrow()) {
        if (subset == SceneSubset::SELECTION) {
            throw std::runtime_error(
                "Reaction selection export is not supported");
        }
        return rdkit_extensions::to_string(*mol_model->getReaction(), format);
    }

    auto mol = extract_mol(mol_model, subset);
    return rdkit_extensions::to_string(*mol, format);
}

void SketcherWidget::addRDKitMolecule(const RDKit::ROMol& mol)
{
    m_mol_model->addMol(mol);
}

void SketcherWidget::addRDKitReaction(const RDKit::ChemicalReaction& rxn)
{
    m_mol_model->addReaction(rxn);
}

boost::shared_ptr<RDKit::ROMol> SketcherWidget::getRDKitMolecule() const
{
    return extract_mol(m_mol_model, SceneSubset::ALL);
}

boost::shared_ptr<RDKit::ChemicalReaction>
SketcherWidget::getRDKitReaction() const
{
    return m_mol_model->getReaction();
}

void SketcherWidget::addFromString(const std::string& text, Format format)
{
    auto probably_a_reaction = [](const auto& text) {
        return boost::starts_with(text, "$RXN") || boost::contains(text, ">>");
    };

    try {
        addRDKitMolecule(*text_to_mol(text, format));
    } catch (const std::exception&) {
        try { // if molecule parsing fails, see if it's a reaction
            addRDKitReaction(*to_rdkit_reaction(text, format));
            return;
        } catch (const std::exception&) {
            // parsing this text as a molecule and as a reaction have both
            // failed, so try to throw the more-relevant exception
            if (probably_a_reaction(text)) {
                throw;
            }
        }
        throw;
    }
}

std::string SketcherWidget::getString(Format format) const
{
    return extract_string(m_mol_model, format, SceneSubset::ALL);
}

void SketcherWidget::cut(Format format)
{
    copy(format, SceneSubset::SELECTION);
    m_mol_model->removeSelected();
}

/**
 * @internal
 * Cut/Copy are currently the only way of extracting a partial molecule from
 * the sketcher widget; the other methods will always extract the entirety of
 * the scene independent of the current selection.
 */
void SketcherWidget::copy(Format format, SceneSubset subset)
{
    std::string text;
    try {
        text = extract_string(m_mol_model, format, subset);
    } catch (const std::exception& exc) {
        show_error_dialog("Copy Error", exc.what(), window());
        return;
    }

    auto data = new QMimeData;
    data->setText(QString::fromStdString(text));
    QApplication::clipboard()->setMimeData(data);

    // SKETCH-2091: Add image content to the clipboard; blocked by SKETCH-1975

#ifdef __EMSCRIPTEN__
    // Use the browser's aync clipboard api to enable copy for the wasm build
    emscripten::val navigator = emscripten::val::global("navigator");
    navigator["clipboard"].call<emscripten::val>("writeText",
                                                 emscripten::val(text));
#endif
}

/**
 * @internal
 * paste is agnostic of NEW_STRUCTURES_REPLACE_CONTENT
 */
void SketcherWidget::paste()
{
    auto data = QApplication::clipboard()->mimeData();
    if (data->hasText()) {
        auto text = data->text().toStdString();
        try {
            addFromString(text, Format::AUTO_DETECT);
        } catch (const std::exception& exc) {
            show_error_dialog("Paste Error", exc.what(), window());
        }
    }
}

void SketcherWidget::importText(const std::string& text, Format format)
{
    if (m_sketcher_model->getValueBool(
            ModelKey::NEW_STRUCTURES_REPLACE_CONTENT)) {
        m_mol_model->clear();
    }

    try {
        addFromString(text, format);
    } catch (const std::exception& exc) {
        show_error_dialog("Import Error", exc.what(), window());
    }
}

void SketcherWidget::showFileExportDialog()
{
    auto dialog = new FileExportDialog(m_sketcher_model, window());
    connect(dialog, &FileExportDialog::exportTextRequested, this,
            [this](Format format) {
                return QString::fromStdString(getString(format));
            });
    dialog->show();
}

void SketcherWidget::showFileSaveImageDialog()
{
    auto dialog = new FileSaveImageDialog(m_sketcher_model, window());
    connect(dialog, &FileSaveImageDialog::exportImageRequested, this,
            [this](auto format, const auto& opts) {
                // SKETCH-1975: this call uses the existing sketcherScene class;
                // update to render directly from this Scene instance
                return get_image_bytes(*m_mol_model->getMol(), format, opts);
            });
    dialog->show();
}

void SketcherWidget::updateWatermark()
{
    bool is_empty = m_sketcher_model->sceneIsEmpty();
    if (is_empty) {
        auto center =
            m_ui->view->mapToScene(m_ui->view->viewport()->rect().center());
        // Center Watermark on view
        m_watermark_item->setPos(
            center.x() - (m_watermark_item->boundingRect().width()) / 2,
            center.y() - (m_watermark_item->boundingRect().height()) / 2);
    }

    m_watermark_item->setVisible(is_empty);
}

void SketcherWidget::onModelValuePinged(ModelKey key, QVariant value)
{
    if (!m_mol_model->hasSelection()) {
        return;
    }
    switch (key) {
        case ModelKey::ELEMENT: {
            auto element = value.value<Element>();
            m_mol_model->mutateSelectedAtoms(element);
            break;
        }
        case ModelKey::ATOM_QUERY: {
            auto query_func = ATOM_TOOL_QUERY_MAP.at(value.value<AtomQuery>());
            auto atom_query =
                std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY>(
                    query_func());
            m_mol_model->mutateSelectedAtoms(atom_query);
            break;
        }
        case ModelKey::BOND_TOOL: {
            auto bond_tool = value.value<BondTool>();
            if (BOND_TOOL_BOND_MAP.count(bond_tool)) {
                auto [bond_type, bond_dir, flippable, cursor_hint_path] =
                    BOND_TOOL_BOND_MAP.at(bond_tool);
                m_mol_model->mutateSelectedBonds(bond_type, bond_dir,
                                                 flippable);
            } else if (BOND_TOOL_QUERY_MAP.count(bond_tool)) {
                auto [query_type, query_func] =
                    BOND_TOOL_QUERY_MAP.at(bond_tool);
                auto bond_query =
                    std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>(
                        query_func());
                bond_query->setTypeLabel(query_type);
                m_mol_model->mutateSelectedBonds(bond_query);
            }
            // Ignore BondTool::ATOM_CHAIN, which should never be pinged
            break;
        }
        case ModelKey::CHARGE_TOOL: {
            auto charge_tool = value.value<ChargeTool>();
            int increment_by = charge_tool == ChargeTool::INCREASE ? 1 : -1;
            m_mol_model->adjustChargeOnSelectedAtoms(increment_by);
            break;
        }
        case ModelKey::DRAW_TOOL: {
            auto tool = value.value<DrawTool>();
            if (tool == DrawTool::ERASE) {
                m_mol_model->removeSelected();
            } else if (tool == DrawTool::EXPLICIT_H) {
                m_mol_model->toggleExplicitHsOnSelectedAtoms();
            }
            break;
        }
        default:
            break;
    }
}

void SketcherWidget::showContextMenu(
    QGraphicsSceneMouseEvent* event,
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups)
{
    QMenu* menu = nullptr;
    if (sgroups.size()) {
        throw std::runtime_error("sgroup context menu not implemented");
    } else if (atoms.size() && bonds.size()) {
        menu = m_selection_context_menu;
    } else if (atoms.size()) {
        menu = m_atom_context_menu;
    } else if (bonds.size()) {
        menu = m_bond_context_menu;
    } else {
        menu = m_background_context_menu;
    }

    menu->move(event->screenPos());
    auto screen_rect = QApplication::screenAt(QCursor::pos())->geometry();
    auto menu_rectangle = menu->geometry();

    // Make sure the menu is not off the screen (or on a different screen)
    if (menu_rectangle.left() < screen_rect.left()) {
        menu_rectangle.moveLeft(screen_rect.left());
    }
    if (menu_rectangle.top() < screen_rect.top()) {
        menu_rectangle.moveTop(screen_rect.top());
    }
    if (menu_rectangle.right() > screen_rect.right()) {
        menu_rectangle.moveRight(screen_rect.right());
    }
    if (menu_rectangle.bottom() > screen_rect.bottom()) {
        menu_rectangle.moveBottom(screen_rect.bottom());
    }
    menu->move(menu_rectangle.topLeft());
    menu->show();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"