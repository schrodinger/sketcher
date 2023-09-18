#include "schrodinger/sketcher/sketcher_widget.h"

#include <QClipboard>
#include <QGraphicsPixmapItem>
#include <QMimeData>
#include <QWidget>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/dialog/file_save_image_dialog.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/scene.h"
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
    m_mol_model(new MolModel(m_undo_stack, this)),
    m_sketcher_model(new SketcherModel(this)),
    m_scene(new Scene(m_mol_model, m_sketcher_model, this))
{
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

    // TODO: Testing components to be removed from the widget
    m_ui->font_sb->setValue(DEFAULT_FONT_SIZE);
    connect(m_ui->font_sb, &QDoubleSpinBox::valueChanged, m_scene,
            &Scene::setFontSize);
    m_ui->carbon_labels_combo->addItems(
        {"No carbon labels", "Terminal carbons only", "All atoms labeled"});
    connect(m_ui->carbon_labels_combo, &QComboBox::currentIndexChanged,
            [this](auto i) { m_scene->setCarbonsLabeled(CarbonLabels(i)); });
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
        m_mol_model->addMolFromText(text, format);
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
                return QString::fromStdString(m_mol_model->getMolText(format));
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

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"