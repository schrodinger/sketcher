#include "schrodinger/sketcher/sketcher_widget.h"

#include <QCursor>
#include <QClipboard>
#include <QGraphicsPixmapItem>
#include <QMimeData>
#include <QWidget>
#include <GraphMol/ROMol.h>
#include <GraphMol/ChemReactions/Reaction.h>
#include <boost/algorithm/string.hpp>

#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/dialog/file_save_image_dialog.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/molviewer/view.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/molops.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
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
    return boost::make_shared<RDKit::ROMol>(*m_mol_model->getMol());
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
    if (m_sketcher_model->hasReaction()) {
        return rdkit_extensions::to_string(*getRDKitReaction(), format);
    }
    return rdkit_extensions::to_string(*getRDKitMolecule(), format);
}

void SketcherWidget::importText(const std::string& text, Format format)
{
    if (m_sketcher_model->getValueBool(
            ModelKey::NEW_STRUCTURES_REPLACE_CONTENT)) {
        m_mol_model->clear();
    }

    auto show_import_failure = [&](const auto& exc) {
        auto text = QString("Import Failed: ") + exc.what();
        show_error_dialog("Import Error", text, window());
    };

    try {
        addFromString(text, format);
    } catch (const std::invalid_argument& exc) {
        show_import_failure(exc);
    } catch (const std::runtime_error& exc) {
        show_import_failure(exc);
    }
}

void SketcherWidget::paste()
{
    auto data = QApplication::clipboard()->mimeData();
    if (data->hasText()) {
        importText(data->text().toStdString(), Format::AUTO_DETECT);
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

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"