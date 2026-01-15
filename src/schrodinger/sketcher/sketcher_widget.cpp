#include "schrodinger/sketcher/sketcher_widget.h"

#include <fmt/format.h>

#include <QApplication>
#include <QClipboard>
#include <QCursor>
#include <QFontDatabase>
#include <QGraphicsPixmapItem>
#include <QKeyEvent>
#include <QMimeData>
#include <QScreen>
#include <QWidget>
#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <rdkit/RDGeneral/Invariant.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/val.h>
#endif

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/dialog/bracket_subgroup_dialog.h"
#include "schrodinger/sketcher/dialog/edit_atom_properties.h"
#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/dialog/file_save_image_dialog.h"
#include "schrodinger/sketcher/dialog/rendering_settings_dialog.h"
#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/menu/atom_context_menu.h"
#include "schrodinger/sketcher/menu/background_context_menu.h"
#include "schrodinger/sketcher/menu/bond_context_menu.h"
#include "schrodinger/sketcher/menu/bracket_subgroup_context_menu.h"
#include "schrodinger/sketcher/menu/selection_context_menu.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/non_molecular_item.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/view.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

/**
 * A simple data class for storing all types of model objects
 */
class ModelObjsByType
{
  public:
    ModelObjsByType(
        const std::unordered_set<const RDKit::Atom*> atoms,
        const std::unordered_set<const RDKit::Bond*> bonds,
        const std::unordered_set<const RDKit::Bond*> secondary_connections,
        const std::unordered_set<const RDKit::SubstanceGroup*> sgroups,
        const std::unordered_set<const NonMolecularObject*>
            non_molecular_objects);
    const std::unordered_set<const RDKit::Atom*> atoms;
    const std::unordered_set<const RDKit::Bond*> bonds;
    const std::unordered_set<const RDKit::Bond*> secondary_connections;
    const std::unordered_set<const RDKit::SubstanceGroup*> sgroups;
    const std::unordered_set<const NonMolecularObject*> non_molecular_objects;
};

ModelObjsByType::ModelObjsByType(
    const std::unordered_set<const RDKit::Atom*> atoms,
    const std::unordered_set<const RDKit::Bond*> bonds,
    const std::unordered_set<const RDKit::Bond*> secondary_connections,
    const std::unordered_set<const RDKit::SubstanceGroup*> sgroups,
    const std::unordered_set<const NonMolecularObject*> non_molecular_objects) :
    atoms(atoms),
    bonds(bonds),
    secondary_connections(secondary_connections),
    sgroups(sgroups),
    non_molecular_objects(non_molecular_objects)
{
}

SketcherWidget::SketcherWidget(QWidget* parent,
                               const InterfaceTypeType interface_type) :
    QWidget(parent),
    m_undo_stack(new QUndoStack(this)),
    m_mol_model(new MolModel(m_undo_stack))
{
    // The tools in ~Scene will access the underlying mol, so we need to
    // make sure the mol model still exists when the scene is destroyed.
    // This is controlled by the order in which parent relationships
    // are defined.
    m_mol_model->setParent(this);

    m_ui.reset(new Ui::SketcherWidgetForm());
    m_ui->setupUi(this);

    // Load fonts BEFORE creating objects that use them
    for (const auto& font_path : {":/resources/fonts/Arimo-Regular.ttf",
                                  ":/resources/fonts/Arimo-Bold.ttf",
                                  ":/resources/fonts/Arimo-Italic.ttf",
                                  ":/resources/fonts/Arimo-BoldItalic.ttf"}) {
        if (QFontDatabase::addApplicationFont(font_path) == -1) {
            throw std::runtime_error(
                fmt::format("Failed to load font: {}", font_path));
        }
    }

    // Now create Scene and other objects that depend on fonts
    m_sketcher_model = new SketcherModel(this);
    m_scene = new Scene(m_mol_model, m_sketcher_model, this);

    m_ui->top_bar_wdg->setModel(m_sketcher_model);
    m_ui->side_bar_wdg->setModel(m_sketcher_model);

    m_ui->view->setScene(m_scene);
    m_ui->view->setMolModel(m_mol_model);
    m_ui->view->setSketcherModel(m_sketcher_model);
    connect(m_scene, &Scene::importTextRequested, this,
            &SketcherWidget::importText);
    connect(m_scene, &Scene::showContextMenuRequested, this,
            &SketcherWidget::showContextMenu);

    // Connect the scene to the view
    connect(m_scene, &Scene::viewportTranslationRequested, m_ui->view,
            &View::translateViewportFromScreenCoords);
    connect(m_scene, &Scene::newCursorHintRequested, m_ui->view,
            &View::onNewCursorHintRequested);
    connect(m_scene, &Scene::atomHovered, this, &SketcherWidget::onAtomHovered);
    connect(m_scene, &Scene::bondHovered, this, &SketcherWidget::onBondHovered);

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

    connect(m_sketcher_model, &SketcherModel::backgroundColorChanged, this,
            &SketcherWidget::onBackgroundColorChanged);

    connect(m_mol_model, &MolModel::selectionChanged, this,
            &SketcherWidget::selectionChanged);
    connect(m_mol_model, &MolModel::modelChanged, this,
            [this](auto what_changed) {
                onMolModelChanged(what_changed & WhatChanged::MOLECULE);
            });
    connect(m_mol_model, &MolModel::coordinatesChanged, [this]() {
        if (!m_ui->view->isDuringPinchGesture() &&
            !m_scene->isDuringAtomDrag()) {
            // if we're in the middle of a mouse drag or a trackpad gesture,
            // wait until that's finished to emit the signal
            emit representationChanged();
        }
    });
    connect(m_ui->view, &View::pinchGestureFinished, this,
            &SketcherWidget::representationChanged);
    connect(m_scene, &Scene::representationChangingAtomDragFinished, this,
            &SketcherWidget::representationChanged);

    connectTopBarSlots();
    connectSideBarSlots();

    // Create and connect the context menus
    m_atom_context_menu =
        new AtomContextMenu(m_sketcher_model, m_mol_model, this);
    m_bond_context_menu = new BondContextMenu(this);
    m_selection_context_menu =
        new SelectionContextMenu(m_sketcher_model, m_mol_model, this);
    m_sgroup_context_menu = new BracketSubgroupContextMenu(this);
    m_background_context_menu =
        new BackgroundContextMenu(m_sketcher_model, this);
    connectContextMenu(*m_atom_context_menu);
    connectContextMenu(*m_bond_context_menu);
    connectContextMenu(*m_selection_context_menu);
    connectContextMenu(*m_sgroup_context_menu);
    connectContextMenu(*m_background_context_menu);

    // create the file and image export dialogs
    m_file_export_dialog = new FileExportDialog(m_sketcher_model, window());
    connect(m_file_export_dialog, &FileExportDialog::exportTextRequested, this,
            [this](Format format) {
                return QString::fromStdString(getString(format));
            });

    m_file_save_image_dialog =
        new FileSaveImageDialog(m_sketcher_model, window());
    connect(m_file_save_image_dialog,
            &FileSaveImageDialog::exportImageRequested, this,
            [this](auto format, const auto& opts) {
                // opts here contains meaningful values only for width and
                // height and background transparency, everything else is
                // default values. We don't call
                // m_sketcher_model->loadRenderOptions because m_sketcher_model
                // already contains the meaningful values (SKETCH-1922)
                return get_image_bytes(*m_scene, format, opts);
            });

    // create the rendering preferences dialog
    m_rendering_settings_dialog =
        new RenderingSettingsDialog(m_sketcher_model, window());

    // force the scene to update the view's cursor now that all of the signals
    // are connected
    m_scene->requestCursorHintUpdate();

    // Update stylesheet and fonts
    setStyleSheet(schrodinger::sketcher::SKETCHER_WIDGET_STYLE);

    // Set up the watermark after loading fonts because the SVG uses them
    m_watermark_item = new QGraphicsPixmapItem();
    m_watermark_item->setPixmap(QPixmap(":icons/2D-Sketcher-watermark.svg"));
    m_watermark_item->setFlag(QGraphicsItem::ItemIgnoresTransformations, true);
    m_scene->addItem(m_watermark_item);
    connect(m_scene, &Scene::changed, this, &SketcherWidget::updateWatermark);
    connect(m_ui->view, &View::resized, this, &SketcherWidget::updateWatermark);

    setInterfaceType(interface_type);
}

SketcherWidget::~SketcherWidget() = default;

/**
 * @internal
 * @param mol_model the model to extract the molecule from
 * @param subset whether to extract everything or just the current selection
 */
static boost::shared_ptr<RDKit::ROMol> extract_mol(MolModel* mol_model,
                                                   SceneSubset subset)
{
    if (subset == SceneSubset::SELECTION) {
        return mol_model->getSelectedMolForExport();
    }
    return mol_model->getMolForExport();
}

/**
 * @internal
 * @param mol_model the model to extract the molecule from
 * @param sketcher_model the sketcher model (for checking if all items selected)
 * @param format the format to serialize the molecule to
 * @param subset whether to extract everything or just the current selection
 */
static std::string extract_string(MolModel* mol_model,
                                  SketcherModel* sketcher_model, Format format,
                                  SceneSubset subset)
{
    if (format == Format::MDL_MOLV2000) {
        throw std::invalid_argument(
            "Sketcher does not support exporting MDL_MOLV2000 format");
    }

    if (mol_model->hasReactionArrow() &&
        (subset == SceneSubset::ALL || sketcher_model->allItemsSelected())) {
        return rdkit_extensions::to_string(*mol_model->getReactionForExport(),
                                           format);
    }

    // Prevent export of partial selection if any non-molecular objects selected
    if (!mol_model->getSelectedNonMolecularObjects().empty()) {
        throw std::runtime_error("Cannot export a partial reaction");
    }

    auto mol = extract_mol(mol_model, subset);
    // remove the 3D conformer if one is present, since we don't want to
    // serialize that
    RDKit::ROMol mol_copy(*mol);
    auto* default_conf = new RDKit::Conformer(mol_copy.getConformer());
    mol_copy.clearConformers();
    mol_copy.addConformer(default_conf);

    return rdkit_extensions::to_string(mol_copy, format);
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
    if (m_copy_of_mol_model_mol == nullptr) {
        m_copy_of_mol_model_mol = m_mol_model->getMolForExport();
    }
    return m_copy_of_mol_model_mol;
}

boost::shared_ptr<RDKit::ChemicalReaction>
SketcherWidget::getRDKitReaction() const
{
    return m_mol_model->getReactionForExport();
}

void SketcherWidget::addFromString(const std::string& text, Format format)
{
    addTextToMolModel(text, format);
}

std::string SketcherWidget::getString(Format format) const
{
    return extract_string(m_mol_model, m_sketcher_model, format,
                          SceneSubset::ALL);
}

QByteArray SketcherWidget::getImageBytes(ImageFormat format) const
{
    RenderOptions renderOptions;
    return get_image_bytes(*m_scene, format, renderOptions);
}

void SketcherWidget::clear()
{
    m_mol_model->clear();
}

bool SketcherWidget::isEmpty() const
{
    return m_mol_model->isEmpty();
}

void SketcherWidget::activateSelectOnlyMode(const SelectionTool tool)
{
    m_sketcher_model->setSelectOnlyModeActive(true);
    setToolbarsVisible(false);
    m_sketcher_model->setSelectToolAllowedWhenSceneEmpty(true);
    m_sketcher_model->setValues(
        {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::SELECT)},
         {ModelKey::SELECTION_TOOL, QVariant::fromValue(tool)}});

    // This avoids unnecessary signals/slots when the palettes are hidden
    // In 2D Overlay, there is a lag when selecting atoms that we
    // fix here by disconnecting signals to updateWidgetsEnabled()
    // slots of each widget in the sketcher.
    m_ui->side_bar_wdg->disconnectAllUpdateWidgetsEnabled();
}

void SketcherWidget::setColorScheme(const ColorScheme color_scheme)
{
    m_sketcher_model->setColorScheme(color_scheme);
}

static std::unordered_set<const RDKit::Atom*>
get_corresponding_atoms_from_different_mol(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const RDKit::ROMol* mol)
{
    std::unordered_set<const RDKit::Atom*> new_atoms;
    std::transform(atoms.begin(), atoms.end(),
                   std::inserter(new_atoms, new_atoms.begin()),
                   [&mol](const auto* atom) {
                       return mol->getAtomWithIdx(atom->getIdx());
                   });
    return new_atoms;
}

static std::unordered_set<const RDKit::Bond*>
get_corresponding_bonds_from_different_mol(
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const RDKit::ROMol* mol)
{
    std::unordered_set<const RDKit::Bond*> new_bonds;
    std::transform(bonds.begin(), bonds.end(),
                   std::inserter(new_bonds, new_bonds.begin()),
                   [&mol](const auto* bond) {
                       return mol->getBondWithIdx(bond->getIdx());
                   });
    return new_bonds;
}

void SketcherWidget::select(const QSet<const RDKit::Atom*>& qatoms,
                            const QSet<const RDKit::Bond*>& qbonds,
                            const SelectMode select_mode)
{
    // convert from qsets to unordered_sets
    std::unordered_set<const RDKit::Atom*> atoms(qatoms.begin(), qatoms.end());
    std::unordered_set<const RDKit::Bond*> bonds(qbonds.begin(), qbonds.end());

    // MolModel expects atoms and bonds that belong to its internal molecule, so
    // we need to get the corresponding atoms and bonds
    const auto* mol = m_mol_model->getMol();
    std::unordered_set<const RDKit::Atom*> mol_model_atoms;
    std::unordered_set<const RDKit::Bond*> mol_model_bonds;
    try {
        mol_model_atoms =
            get_corresponding_atoms_from_different_mol(atoms, mol);
        mol_model_bonds =
            get_corresponding_bonds_from_different_mol(bonds, mol);
    } catch (const Invar::Invariant&) {
        throw std::runtime_error(
            "Atoms and bonds must belong to the Sketcher molecule.");
    }
    m_mol_model->select(mol_model_atoms, mol_model_bonds, {}, {}, {},
                        select_mode);
}

void SketcherWidget::clearSelection()
{
    m_mol_model->clearSelection();
}

void SketcherWidget::selectAll()
{
    m_mol_model->selectAll();
}

QSet<const RDKit::Atom*> SketcherWidget::getSelectedAtoms() const
{
    auto atoms = m_mol_model->getSelectedAtoms();
    // Python won't respect the constness of the Atom pointers, so we return
    // atoms from a copy of the molecule instead
    auto mol = getRDKitMolecule();
    auto atoms_from_copy =
        get_corresponding_atoms_from_different_mol(atoms, mol.get());
    return QSet(atoms_from_copy.begin(), atoms_from_copy.end());
}

QSet<const RDKit::Bond*> SketcherWidget::getSelectedBonds() const
{
    auto bonds = m_mol_model->getSelectedBonds();
    // Python won't respect the constness of the Bond pointers, so we return
    // bonds from a copy of the molecule instead
    auto mol = getRDKitMolecule();
    auto bonds_from_copy =
        get_corresponding_bonds_from_different_mol(bonds, mol.get());
    return QSet(bonds_from_copy.begin(), bonds_from_copy.end());
}

void SketcherWidget::fitToScreen(bool selection_only)
{
    m_ui->view->fitToScreen(selection_only);
}

void SketcherWidget::setInterfaceType(InterfaceTypeType interface_type)
{
    m_sketcher_model->setValue(ModelKey::INTERFACE_TYPE, interface_type);
}

std::string SketcherWidget::getClipboardContents() const
{
    auto data = QApplication::clipboard()->mimeData();
    if (data->hasText()) {
        return data->text().toStdString();
    }
    return "";
}

void SketcherWidget::setClipboardContents(std::string text) const
{
    auto data = new QMimeData;
    data->setText(QString::fromStdString(text));
    QApplication::clipboard()->setMimeData(data);
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
        text = extract_string(m_mol_model, m_sketcher_model, format, subset);
    } catch (const std::exception& exc) {
        show_error_dialog("Copy Error", exc.what(), window());
        return;
    }
    setClipboardContents(text);
    // SKETCH-2091: Add image content to the clipboard; blocked by SKETCH-1975

#ifdef __EMSCRIPTEN__
    // Use the browser's aync clipboard api to enable copy for the wasm build
    emscripten::val navigator = emscripten::val::global("navigator");
    navigator["clipboard"].call<emscripten::val>("writeText",
                                                 emscripten::val(text));
#endif
}

void SketcherWidget::copyAsImage()
{
    QByteArray image_bytes;
    RenderOptions opts;
    image_bytes = get_image_bytes(*m_scene, ImageFormat::PNG, opts);
    QImage image;
    image.loadFromData(image_bytes, "PNG");
    QApplication::clipboard()->setImage(image);
}

/**
 * @internal
 * paste is agnostic of NEW_STRUCTURES_REPLACE_CONTENT
 */
void SketcherWidget::pasteAt(std::optional<QPointF> position)
{
    auto text = getClipboardContents();
    if (text.empty()) {
        return;
    }
    // On WASM builds, RDKit doesn't like Windows newline characters, so we
    // explicitly remove the /r's, which converts Windows-style newlines to
    // Unix-style
    std::erase_if(text, [](char c) { return c == '\r'; });
    std::optional<RDGeom::Point3D> mol_position = std::nullopt;
    if (position.has_value()) {
        // Convert the position from global to scene coordinates
        auto scene_position = m_ui->view->mapToScene(
            m_ui->view->mapFromGlobal(position.value()).toPoint());
        // and then to mol coordinates
        mol_position = to_mol_xy(scene_position);
    }
    try {
        addTextToMolModel(text, Format::AUTO_DETECT, mol_position,
                          /*recenter_view*/ false);
    } catch (const std::exception& exc) {
        show_error_dialog("Paste Error", exc.what(), window());
    }
}

/**
 * @internal
 * paste is agnostic of NEW_STRUCTURES_REPLACE_CONTENT
 */
void SketcherWidget::paste()
{
    pasteAt(std::nullopt);
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
    m_file_export_dialog->show();
}

void SketcherWidget::showFileSaveImageDialog()
{
    m_file_save_image_dialog->show();
}

void SketcherWidget::showRenderingSettingsDialog()
{
    m_rendering_settings_dialog->show();
}

void SketcherWidget::showEditAtomPropertiesDialog(
    const RDKit::Atom* const atom, const bool set_to_allowed_list)
{
    auto dialog = new EditAtomPropertiesDialog(atom, m_mol_model, this);
    if (set_to_allowed_list) {
        dialog->switchToAllowedList();
    }
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

void SketcherWidget::showBracketSubgroupDialogForAtoms(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    auto dialog = new BracketSubgroupDialog(m_mol_model, this);
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setAtoms(atoms);
    dialog->show();
}

void SketcherWidget::showBracketSubgroupDialogForSGroup(
    const RDKit::SubstanceGroup* const s_group)
{
    auto dialog = new BracketSubgroupDialog(m_mol_model, this);
    dialog->setAttribute(Qt::WA_DeleteOnClose);
    dialog->setSubgroup(s_group);
    dialog->show();
}

void SketcherWidget::connectTopBarSlots()
{
    connect(m_ui->top_bar_wdg, &SketcherTopBar::undoRequested, m_undo_stack,
            &QUndoStack::undo);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::redoRequested, m_undo_stack,
            &QUndoStack::redo);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::cleanupRequested,
            [this](bool selection_only) {
                if (selection_only) {
                    m_mol_model->cleanUpSelection();
                } else {
                    m_mol_model->regenerateCoordinates();
                }
                m_ui->view->fitToScreen(selection_only);
            });
    connect(m_ui->top_bar_wdg, &SketcherTopBar::fitToScreenRequested,
            m_ui->view, &View::fitToScreen);

    // Connect "More Actions" menu
    connect(m_ui->top_bar_wdg, &SketcherTopBar::flipHorizontalRequested,
            m_mol_model, &MolModel::flipAllHorizontal);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::flipVerticalRequested,
            m_mol_model, &MolModel::flipAllVertical);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::aromatizeRequested, m_mol_model,
            &MolModel::aromatize);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::kekulizeRequested, m_mol_model,
            &MolModel::kekulize);
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
    connect(m_ui->top_bar_wdg,
            &SketcherTopBar::adjustRenderingSettingsRequested, this,
            &SketcherWidget::showRenderingSettingsDialog);
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
    using RDKitAtoms = std::unordered_set<const RDKit::Atom*>;
    connect(
        &menu, &ModifyAtomsMenu::requestElementChange, m_mol_model,
        qOverload<const RDKitAtoms&, const Element>(&MolModel::mutateAtoms));
    connect(&menu, &ModifyAtomsMenu::adjustChargeRequested, m_mol_model,
            &MolModel::adjustChargeOnAtoms);

    connect(&menu, &ModifyAtomsMenu::addRemoveExplicitHydrogensRequested,
            m_mol_model, &MolModel::toggleExplicitHsOnAtoms);
    connect(&menu, &ModifyAtomsMenu::adjustRadicalElectronsRequested,
            m_mol_model, &MolModel::adjustRadicalElectronsOnAtoms);

    connect(&menu, &ModifyAtomsMenu::showEditAtomPropertiesRequested, this,
            &SketcherWidget::showEditAtomPropertiesDialog);
    connect(
        &menu, &ModifyAtomsMenu::changeTypeRequested, m_mol_model,
        qOverload<const RDKitAtoms&, const AtomQuery>(&MolModel::mutateAtoms));
    connect(&menu, &ModifyAtomsMenu::newRGroupRequested, m_mol_model,
            qOverload<const RDKitAtoms&>(&MolModel::mutateRGroups));
    connect(&menu, &ModifyAtomsMenu::existingRGroupRequested, m_mol_model,
            qOverload<const RDKitAtoms&, const unsigned int>(
                &MolModel::mutateRGroups));

    if (auto context_menu = dynamic_cast<const AtomContextMenu*>(&menu)) {
        connect(context_menu, &AtomContextMenu::bracketSubgroupDialogRequested,
                this, &SketcherWidget::showBracketSubgroupDialogForAtoms);
        connect(
            context_menu, &AtomContextMenu::deleteRequested, this,
            [this](auto atoms) { m_mol_model->remove(atoms, {}, {}, {}, {}); });
    }
}

void SketcherWidget::connectContextMenu(const ModifyBondsMenu& menu)
{
    connect(&menu, &ModifyBondsMenu::changeTypeRequested, this,
            [this](auto bond_tool, auto bonds) {
                m_mol_model->mutateBonds(bonds, bond_tool);
            });
    connect(&menu, &ModifyBondsMenu::changeQueryRequested, this,
            [this](auto bond_tool, auto bonds) {
                m_mol_model->setBondTopology(bonds, bond_tool);
            });
    connect(&menu, &ModifyBondsMenu::flipRequested, this, [this](auto bonds) {
        // flip substituent only makes sense if there is only one bond
        if (bonds.size() != 1) {
            throw std::runtime_error(
                "Cannot flip substituent for multiple bonds");
        }
        m_mol_model->flipSubstituent(*(bonds.begin()));
    });
    if (auto context_menu = dynamic_cast<const BondContextMenu*>(&menu)) {
        connect(context_menu, &BondContextMenu::deleteRequested, this,
                [this](auto bonds, auto secondary_conditions) {
                    m_mol_model->remove({}, bonds, secondary_conditions, {},
                                        {});
                });
    }
}
void SketcherWidget::connectContextMenu(const SelectionContextMenu& menu)
{
    connect(&menu, &SelectionContextMenu::cleanUpRegionRequested, m_mol_model,
            &MolModel::cleanUpSelection);
    connect(&menu, &SelectionContextMenu::invertSelectionRequested, m_mol_model,
            &MolModel::invertSelection);
    connect(&menu, &SelectionContextMenu::cutRequested, this,
            &SketcherWidget::cut);
    connect(&menu, &SelectionContextMenu::copyRequested, this,
            &SketcherWidget::copy);
    connect(&menu, &SelectionContextMenu::copyAsImageRequested, this,
            &SketcherWidget::copyAsImage);
    connect(&menu, &SelectionContextMenu::flipRequested, m_mol_model,
            &MolModel::flipSelection);
    connect(&menu, &SelectionContextMenu::flipHorizontalRequested, m_mol_model,
            &MolModel::flipSelectionHorizontal);
    connect(&menu, &SelectionContextMenu::flipVerticalRequested, m_mol_model,
            &MolModel::flipSelectionVertical);
    connectContextMenu(*menu.m_modify_atoms_menu);
    connectContextMenu(*menu.m_modify_bonds_menu);
    connect(&menu, &SelectionContextMenu::bracketSubgroupDialogRequested, this,
            &SketcherWidget::showBracketSubgroupDialogForAtoms);
    connect(&menu, &SelectionContextMenu::variableAttachmentBondRequested, this,
            [this](const auto& atoms) {
                auto undo_raii = m_mol_model->createUndoMacro(
                    "Add variable attachment bond");
                m_mol_model->addVariableAttachmentBond(atoms);
                m_mol_model->clearSelection();
            });
    connect(&menu, &SelectionContextMenu::deleteRequested, this,
            [this](auto atoms, auto bonds, auto secondary_connections,
                   auto sgroups, auto non_mol_objs) {
                m_mol_model->remove(atoms, bonds, secondary_connections,
                                    sgroups, non_mol_objs);
            });
}

void SketcherWidget::connectContextMenu(const BracketSubgroupContextMenu& menu)
{
    connect(&menu, &BracketSubgroupContextMenu::bracketSubgroupDialogRequested,
            this, [this](auto sgroups) {
                // modification dialog only makes sense for one sgroup
                if (sgroups.size() != 1) {
                    throw std::runtime_error(
                        "Cannot modify more multiple bracket groups");
                }
                showBracketSubgroupDialogForSGroup(*sgroups.begin());
            });
    connect(
        &menu, &BracketSubgroupContextMenu::deleteRequested, this,
        [this](auto sgroups) { m_mol_model->remove({}, {}, {}, sgroups, {}); });
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
    connect(&menu, &BackgroundContextMenu::copyAsImageRequested, this,
            &SketcherWidget::copyAsImage);
    connect(&menu, &BackgroundContextMenu::pasteRequested, this,
            &SketcherWidget::pasteAt);
    connect(&menu, &BackgroundContextMenu::clearRequested, m_mol_model,
            &MolModel::clear);
}

void SketcherWidget::showContextMenu(
    QGraphicsSceneMouseEvent* event,
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const std::unordered_set<const RDKit::Bond*>& secondary_connections,
    const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{
    if (m_sketcher_model && m_sketcher_model->isSelectOnlyModeActive()) {
        // context menus are disabled when select-only mode is active to prevent
        // the user from mutating the structure via the context menu
        return;
    }

    AbstractContextMenu* menu = nullptr;
    bool bond_selected = bonds.size() || secondary_connections.size();
    if (sgroups.size()) {
        menu = m_sgroup_context_menu;
    } else if (atoms.size() && bond_selected) {
        menu = m_selection_context_menu;
    } else if (atoms.size()) {
        menu = m_atom_context_menu;
    } else if (bond_selected) {
        menu = m_bond_context_menu;
    } else if (non_molecular_objects.size()) {
        /** a non molecular object is clicked, do nothing as there's no context
         * menu for it.
         * */
        return;
    } else {
        // show the background context menu
        menu = m_background_context_menu;
    }
    menu->setContextItems(atoms, bonds, secondary_connections, sgroups,
                          non_molecular_objects);

    menu->move(event->screenPos());
    auto screen_rect = QApplication::screenAt(QCursor::pos())->geometry();

    // Ensure the size is updated before checking geometry. This is necessary
    // the first time we execute this because the context menu hasn't been shown
    // yet
    menu->adjustSize();
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

void SketcherWidget::setToolbarsVisible(const bool visible)
{
    m_ui->side_bar_wdg->setVisible(visible);
    m_ui->top_bar_wdg->setVisible(visible);
    // also hide the line that's between the workspace and the top toolbar
    m_ui->line->setVisible(visible);
}

void SketcherWidget::keyPressEvent(QKeyEvent* event)
{
    QWidget::keyPressEvent(event);

    if (m_sketcher_model && m_sketcher_model->isSelectOnlyModeActive()) {
        // keyboard shortcuts are disabled when select-only mode is active to
        // prevent the user from switching tools
        return;
    }

    auto cursor_pos =
        m_ui->view->mapToScene(m_ui->view->mapFromGlobal(QCursor::pos()));
    auto [atoms, bonds, secondary_connections, sgroups, non_molecular_objects] =
        m_scene->getModelObjects(SceneSubset::SELECTED_OR_HOVERED, &cursor_pos);
    ModelObjsByType targets(atoms, bonds, secondary_connections, sgroups,
                            non_molecular_objects);

    bool handled = handleCommonKeyboardShortcuts(event, cursor_pos, targets);
    if (!handled) {
        if (m_sketcher_model->getToolSet() == ToolSet::ATOMISTIC) {
            handleAtomisticKeyboardShortcuts(event, cursor_pos, targets);
        } else if (m_sketcher_model->getMonomerToolType() ==
                   MonomerToolType::AMINO_ACID) {
            handleAminoAcidKeyboardShortcuts(event, cursor_pos, targets);
        } else {
            handleNucleicAcidKeyboardShortcuts(event, cursor_pos, targets);
        }
    }
}

void SketcherWidget::updateModelForKeyboardShortcut(
    const bool has_targets, const std::pair<ModelKey, QVariant>& kv_pair,
    std::unordered_map<ModelKey, QVariant> kv_pairs,
    const ModelObjsByType& targets)
{
    // unless we are working on a selection, switch model tool
    if (!m_mol_model->hasSelection()) {
        kv_pairs.emplace(kv_pair.first, kv_pair.second);
        m_sketcher_model->setValues(kv_pairs);
    }
    // Finally, interact with models
    if (has_targets) {
        applyModelValuePingToTargets(
            kv_pair.first, kv_pair.second, targets.atoms, targets.bonds,
            targets.secondary_connections, targets.sgroups,
            targets.non_molecular_objects);
    }
}

bool SketcherWidget::handleCommonKeyboardShortcuts(
    QKeyEvent* event, const QPointF& cursor_pos, const ModelObjsByType& targets)
{
    switch (Qt::Key(event->key())) {
        case Qt::Key_Backspace:
        case Qt::Key_Delete: {
            std::pair<ModelKey, QVariant> kv_pair = {
                ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::ERASE)};
            bool has_target_atoms = !targets.atoms.empty();
            bool has_target_bonds = !targets.bonds.empty();
            bool has_target_secondary_connections =
                !targets.secondary_connections.empty();
            bool has_target_non_molecular_objects =
                !targets.non_molecular_objects.empty();
            bool has_targets = has_target_atoms || has_target_bonds ||
                               has_target_secondary_connections ||
                               has_target_non_molecular_objects;
            updateModelForKeyboardShortcut(has_targets, kv_pair, {}, targets);
            return true;
        }
        case Qt::Key_Space:
            // Switch back to the last used select tool
            if (!m_sketcher_model->sceneIsEmpty()) {
                m_sketcher_model->setValue(ModelKey::DRAW_TOOL,
                                           DrawTool::SELECT);
            }
            return true;
        default:
            return false;
    }
}

void SketcherWidget::handleAtomisticKeyboardShortcuts(
    QKeyEvent* event, const QPointF& cursor_pos, const ModelObjsByType& targets)
{
    std::pair<ModelKey, QVariant> kv_pair;
    std::unordered_map<ModelKey, QVariant> kv_pairs;
    bool has_targets;
    bool has_target_atoms = !targets.atoms.empty();
    bool has_target_bonds = !targets.bonds.empty();
    // we don't need to worry about secondary_connections here since they don't
    // exist in atomistic models

    switch (Qt::Key(event->key())) {
        case Qt::Key_Minus:
        case Qt::Key_Plus:
        case Qt::Key_Equal: {
            auto tool = event->key() == Qt::Key_Minus ? ChargeTool::DECREASE
                                                      : ChargeTool::INCREASE;
            kv_pair = {ModelKey::CHARGE_TOOL, QVariant::fromValue(tool)};
            kv_pairs = {
                {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::CHARGE)}};
            has_targets = has_target_atoms;
            updateModelForKeyboardShortcut(has_targets, kv_pair, kv_pairs,
                                           targets);
            return;
        }
        case Qt::Key_0:
        case Qt::Key_1:
        case Qt::Key_2:
        case Qt::Key_3: {
            std::unordered_map<int, BondTool> key_to_bond_tool = {
                {Qt::Key_0, BondTool::ZERO},
                {Qt::Key_1, BondTool::SINGLE},
                {Qt::Key_2, BondTool::DOUBLE},
                {Qt::Key_3, BondTool::TRIPLE},
            };
            auto tool = key_to_bond_tool[event->key()];
            kv_pair = {ModelKey::BOND_TOOL, QVariant::fromValue(tool)};
            kv_pairs = {
                {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::BOND)}};
            has_targets = has_target_bonds;
            updateModelForKeyboardShortcut(has_targets, kv_pair, kv_pairs,
                                           targets);
            return;
        }

        case Qt::Key_D: {
            auto to_atom = RDKit::Atom("H");
            to_atom.setIsotope(2);
            m_mol_model->mutateAtoms(targets.atoms, to_atom);
            return;
        }
        case Qt::Key_T: {
            auto to_atom = RDKit::Atom("H");
            to_atom.setIsotope(3);
            m_mol_model->mutateAtoms(targets.atoms, to_atom);
            return;
        }
        /* case Qt::Key_Escape:
            // TODO: SKETCH-1184 SKETCH-2045
        */
        default: {
            // Check whether pressed key corresponds to an element
            auto key_text = event->text().toStdString();
            int atomic_number = -1;
            try {
                atomic_number = symbol_to_atomic_number(key_text);
            } catch (Invar::Invariant const&) {
                // Persist invalid atomic number if not an element
            }
            if (atomic_number != -1 && event->modifiers() == Qt::NoModifier) {
                auto element = static_cast<Element>(atomic_number);
                has_targets = has_target_atoms;
                kv_pair = {ModelKey::ELEMENT, QVariant::fromValue(element)};
                kv_pairs = {
                    {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::ATOM)},
                    {ModelKey::ATOM_TOOL,
                     QVariant::fromValue(AtomTool::ELEMENT)}};
                updateModelForKeyboardShortcut(has_targets, kv_pair, kv_pairs,
                                               targets);
            }
            return;
        }
    }
}

void SketcherWidget::handleAminoAcidKeyboardShortcuts(
    QKeyEvent* event, const QPointF& cursor_pos, const ModelObjsByType& targets)
{
    static const std::unordered_map<Qt::Key, AminoAcidTool> KEY_TO_AMINO_ACID{
        // clang-format off
        {Qt::Key_A, AminoAcidTool::ALA},
        {Qt::Key_R, AminoAcidTool::ARG},
        {Qt::Key_N, AminoAcidTool::ASN},
        {Qt::Key_D, AminoAcidTool::ASP},
        {Qt::Key_C, AminoAcidTool::CYS},
        {Qt::Key_Q, AminoAcidTool::GLN},
        {Qt::Key_E, AminoAcidTool::GLU},
        {Qt::Key_G, AminoAcidTool::GLY},
        {Qt::Key_H, AminoAcidTool::HIS},
        {Qt::Key_I, AminoAcidTool::ILE},
        {Qt::Key_L, AminoAcidTool::LEU},
        {Qt::Key_K, AminoAcidTool::LYS},
        {Qt::Key_M, AminoAcidTool::MET},
        {Qt::Key_F, AminoAcidTool::PHE},
        {Qt::Key_P, AminoAcidTool::PRO},
        {Qt::Key_S, AminoAcidTool::SER},
        {Qt::Key_T, AminoAcidTool::THR},
        {Qt::Key_W, AminoAcidTool::TRP},
        {Qt::Key_Y, AminoAcidTool::TYR},
        {Qt::Key_V, AminoAcidTool::VAL},
        {Qt::Key_X, AminoAcidTool::UNK},
    }; // clang-format on
    auto key = Qt::Key(event->key());
    if (KEY_TO_AMINO_ACID.contains(key)) {
        auto amino_acid_tool = KEY_TO_AMINO_ACID.at(key);
        bool has_targets = !targets.atoms.empty();
        std::pair<ModelKey, QVariant> kv_pair = {
            ModelKey::AMINO_ACID_TOOL, QVariant::fromValue(amino_acid_tool)};
        std::unordered_map<ModelKey, QVariant> kv_pairs = {
            {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
            {ModelKey::MONOMER_TOOL_TYPE,
             QVariant::fromValue(MonomerToolType::AMINO_ACID)},
        };
        updateModelForKeyboardShortcut(has_targets, kv_pair, kv_pairs, targets);
    }
}

void SketcherWidget::handleNucleicAcidKeyboardShortcuts(
    QKeyEvent* event, const QPointF& cursor_pos, const ModelObjsByType& targets)
{
    // behavior for the keyboard keys that represent nucleobases, depending on
    // the currently active tool
    //  - the base name to switch the custom nucleotide tool to
    //  - the base to switch the RNA/DNA nucleotide tool to
    //  - the monomer tool to switch to
    static const std::unordered_map<
        Qt::Key, std::tuple<QString, StdNucleobase, NucleicAcidTool>>
        BASE_KEYS{
            {Qt::Key_A, {"A", StdNucleobase::A, NucleicAcidTool::A}},
            {Qt::Key_C, {"C", StdNucleobase::C, NucleicAcidTool::C}},
            {Qt::Key_G, {"G", StdNucleobase::G, NucleicAcidTool::G}},
            {Qt::Key_U, {"U", StdNucleobase::U_OR_T, NucleicAcidTool::U}},
            {Qt::Key_T, {"T", StdNucleobase::U_OR_T, NucleicAcidTool::T}},
            {Qt::Key_N, {"N", StdNucleobase::N, NucleicAcidTool::N}},
        };

    // behavior for the keyboard keys that represent sugars, depending on
    // the currently active tool
    //  - the sugar name to switch the custom nucleotide tool to
    //  - the RNA/DNA nucleotide tool to switch to
    //  - the monomer tool to switch to
    static const std::unordered_map<
        Qt::Key, std::tuple<QString, NucleicAcidTool, NucleicAcidTool>>
        SUGAR_KEYS{
            {Qt::Key_R,
             {"R", NucleicAcidTool::RNA_NUCLEOTIDE, NucleicAcidTool::R}},
            {Qt::Key_D,
             {"dR", NucleicAcidTool::DNA_NUCLEOTIDE, NucleicAcidTool::dR}},
        };

    // behavior for the P key (i.e. phosphate), depending on
    // the currently active tool
    //  - the phosphate name to switch the custom nucleotide tool to (in case
    //    the user had previously switched to a modified phosphate group)
    //  - the monomer tool to switch to
    // Note that this key has no effect if the RNA/DNA nucleotide tool is
    // active, since those tools don't allow changes to the phosphate group
    static const std::tuple<QString, NucleicAcidTool> P_KEY{"P",
                                                            NucleicAcidTool::P};

    auto key = Qt::Key(event->key());
    std::optional<std::pair<ModelKey, QVariant>> kv_pair;
    auto current_tool = m_sketcher_model->getNucleicAcidTool();
    switch (current_tool) {
        case NucleicAcidTool::RNA_NUCLEOTIDE:
        case NucleicAcidTool::DNA_NUCLEOTIDE:
            if (BASE_KEYS.contains(key)) {
                // switch the base of the current full nucleotide tool
                auto base = std::get<StdNucleobase>(BASE_KEYS.at(key));
                auto model_key = current_tool == NucleicAcidTool::RNA_NUCLEOTIDE
                                     ? ModelKey::RNA_NUCLEOBASE
                                     : ModelKey::DNA_NUCLEOBASE;
                kv_pair = {model_key, QVariant::fromValue(base)};
            } else if (SUGAR_KEYS.contains(key)) {
                // switch the specified full nucleotide tool
                auto tool = std::get<1>(SUGAR_KEYS.at(key));
                kv_pair = {ModelKey::NUCLEIC_ACID_TOOL,
                           QVariant::fromValue(tool)};
            }
            // the P key has no effect in this scenario
            break;
        case NucleicAcidTool::CUSTOM_NUCLEOTIDE: {
            auto nt = m_sketcher_model->getCustomNucleotide();
            if (BASE_KEYS.contains(key)) {
                // switch the base of the current tool
                auto base = std::get<QString>(BASE_KEYS.at(key));
                std::get<1>(nt) = base;
            } else if (SUGAR_KEYS.contains(key)) {
                // switch the sugar of the current tool
                auto sugar = std::get<QString>(SUGAR_KEYS.at(key));
                std::get<0>(nt) = sugar;
            } else if (key == Qt::Key_P) {
                // remove any phosphate modifications from the current tool
                auto phosphate = std::get<QString>(P_KEY);
                std::get<2>(nt) = phosphate;
            } else {
                break;
            }
            kv_pair = {ModelKey::CUSTOM_NUCLEOTIDE, QVariant::fromValue(nt)};
            break;
        }
        default: {
            // a monomer tool is active, so switch to the requested monomer tool
            // (i.e. not a full nucleotide tool)
            NucleicAcidTool tool;
            if (BASE_KEYS.contains(key)) {
                tool = std::get<NucleicAcidTool>(BASE_KEYS.at(key));
            } else if (SUGAR_KEYS.contains(key)) {
                tool = std::get<2>(SUGAR_KEYS.at(key));
            } else if (key == Qt::Key_P) {
                tool = std::get<NucleicAcidTool>(P_KEY);
            } else {
                break;
            }
            kv_pair = {ModelKey::NUCLEIC_ACID_TOOL, QVariant::fromValue(tool)};
        }
    }

    if (kv_pair.has_value()) {
        std::unordered_map<ModelKey, QVariant> kv_pairs = {
            {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MONOMER)},
            {ModelKey::MONOMER_TOOL_TYPE,
             QVariant::fromValue(MonomerToolType::NUCLEIC_ACID)},
        };
        bool has_targets = !targets.atoms.empty();
        updateModelForKeyboardShortcut(has_targets, *kv_pair, kv_pairs,
                                       targets);
    }
}

void SketcherWidget::onBackgroundColorChanged(const QColor& color)
{
    auto view = m_ui->view;
    if (!view) {
        return;
    }
    view->setBackgroundBrush(QBrush(color));
}

void SketcherWidget::onModelValuePinged(ModelKey key, QVariant value)
{
    auto view = m_ui->view;
    if (!view) {
        return;
    }
    auto [atoms, bonds, secondary_connections, sgroups, non_molecular_objects] =
        m_scene->getModelObjects(SceneSubset::SELECTION);
    applyModelValuePingToTargets(key, value, atoms, bonds,
                                 secondary_connections, sgroups,
                                 non_molecular_objects);
}

void SketcherWidget::applyModelValuePingToTargets(
    const ModelKey key, const QVariant value,
    const std::unordered_set<const RDKit::Atom*> atoms,
    const std::unordered_set<const RDKit::Bond*> bonds,
    const std::unordered_set<const RDKit::Bond*> secondary_connections,
    const std::unordered_set<const RDKit::SubstanceGroup*> sgroups,
    const std::unordered_set<const NonMolecularObject*> non_molecular_objects)
{
    switch (key) {
        case ModelKey::ELEMENT: {
            auto element = value.value<Element>();
            m_mol_model->mutateAtoms(atoms, element);
            break;
        }
        case ModelKey::ATOM_QUERY: {
            auto atom_query = value.value<AtomQuery>();
            m_mol_model->mutateAtoms(atoms, atom_query);
            break;
        }
        case ModelKey::BOND_TOOL: {
            auto bond_tool = value.value<BondTool>();
            m_mol_model->mutateBonds(bonds, bond_tool);
            break;
        }
        case ModelKey::CHARGE_TOOL: {
            auto charge_tool = value.value<ChargeTool>();
            int increment_by = charge_tool == ChargeTool::INCREASE ? 1 : -1;
            m_mol_model->adjustChargeOnAtoms(atoms, increment_by);
            break;
        }
        case ModelKey::DRAW_TOOL: {
            auto tool = value.value<DrawTool>();
            if (tool == DrawTool::ERASE) {
                m_mol_model->remove(atoms, bonds, secondary_connections,
                                    sgroups, non_molecular_objects);
            } else if (tool == DrawTool::EXPLICIT_H) {
                m_mol_model->toggleExplicitHsOnAtoms(atoms);
            }
            break;
        }
        default:
            break;
    }
}

void SketcherWidget::onAtomHovered(const RDKit::Atom* atom)
{
    if (atom != nullptr) {
        atom = getRDKitMolecule()->getAtomWithIdx(atom->getIdx());
    }
    emit atomHovered(atom);
}

void SketcherWidget::onBondHovered(const RDKit::Bond* bond)
{
    if (bond != nullptr) {
        bond = getRDKitMolecule()->getBondWithIdx(bond->getIdx());
    }
    emit bondHovered(bond);
}

void SketcherWidget::onMolModelChanged(const bool molecule_changed)
{
    // even if what_changed doesn't contain MOLECULE, it's possible that the
    // molecule's coordinates have changed so we should clear our cached
    // copy of the molecule no matter what
    m_copy_of_mol_model_mol = nullptr;

    if (molecule_changed) {
        // update MOLECULE_TYPE (i.e. is the Sketcher workspace empty,
        // atomistic, or monomeric) and TOOL_SET (i.e. does the side bar
        // currently display atomistic or monomeric tools). Note that we don't
        // allow the workspace to contain both atomistic and monomeric models at
        // the same time; it must be one or the other (or empty).
        std::unordered_map<ModelKey, QVariant> kv_pairs;
        auto interface = m_sketcher_model->getInterfaceType();
        if (!m_mol_model->hasMolecularObjects()) {
            kv_pairs.emplace(ModelKey::MOLECULE_TYPE,
                             QVariant::fromValue(MoleculeType::EMPTY));
        } else if (m_mol_model->isMonomeric()) {
            kv_pairs.emplace(ModelKey::MOLECULE_TYPE,
                             QVariant::fromValue(MoleculeType::MONOMERIC));
            if (interface & InterfaceType::MONOMERIC) {
                kv_pairs.emplace(ModelKey::TOOL_SET,
                                 QVariant::fromValue(ToolSet::MONOMERIC));
            }
        } else {
            kv_pairs.emplace(ModelKey::MOLECULE_TYPE,
                             QVariant::fromValue(MoleculeType::ATOMISTIC));
            if (interface & InterfaceType::ATOMISTIC) {
                kv_pairs.emplace(ModelKey::TOOL_SET,
                                 QVariant::fromValue(ToolSet::ATOMISTIC));
            }
        }
        m_sketcher_model->setValues(kv_pairs);
    }

    if (molecule_changed) {
        emit moleculeChanged();
    } else {
        emit representationChanged();
    }
}

bool SketcherWidget::handleShortcutAction(const QKeySequence& key)
{
    return m_ui->top_bar_wdg->handleShortcutAction(key);
}

QRectF SketcherWidget::getSceneRect()
{
    return m_scene->getSceneItemsBoundingRect();
}

void SketcherWidget::addTextToMolModel(
    const std::string& text, const rdkit_extensions::Format format,
    const std::optional<RDGeom::Point3D> position, const bool recenter_view)
{
    auto mol_or_reaction = convert_text_to_mol_or_reaction(text, format);
    auto* mol_ptr_ptr =
        std::get_if<boost::shared_ptr<RDKit::RWMol>>(&mol_or_reaction);
    // we assume that reactions are atomistic
    bool mol_to_add_is_monomeric =
        (mol_ptr_ptr && rdkit_extensions::isMonomeric(**mol_ptr_ptr));

    auto interface_type = m_sketcher_model->getInterfaceType();
    auto cur_molecule_type = m_sketcher_model->getMoleculeType();
    if (mol_to_add_is_monomeric &&
        !(interface_type & InterfaceType::MONOMERIC)) {
        throw std::runtime_error("Monomeric models not allowed");
    } else if (!mol_to_add_is_monomeric &&
               !(interface_type & InterfaceType::ATOMISTIC)) {
        throw std::runtime_error("Atomistic models not allowed");
    } else if (cur_molecule_type == MoleculeType::ATOMISTIC &&
               mol_to_add_is_monomeric) {
        throw std::runtime_error(
            "Cannot add a monomeric model when the Sketcher already contains "
            "an atomistic model");
    } else if (cur_molecule_type == MoleculeType::MONOMERIC &&
               !mol_to_add_is_monomeric) {
        throw std::runtime_error(
            "Cannot add an atomistic model when the Sketcher already contains "
            "a monomeric model");
    }
    add_mol_or_reaction_to_mol_model(*m_mol_model, mol_or_reaction, position,
                                     recenter_view);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"
