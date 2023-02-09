#include "schrodinger/sketcher/sketcher_widget.h"

#include <QGraphicsPixmapItem>
#include <QWidget>

#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_sketcher_widget.h"

namespace schrodinger
{
namespace sketcher
{

SketcherWidget::SketcherWidget(QWidget* parent) : QWidget(parent)
{
    m_ui.reset(new Ui::SketcherWidgetForm());
    m_ui->setupUi(this);

    m_scene = new Scene(this);
    m_ui->view->setScene(m_scene);

    m_sketcher_model = new SketcherModel(this);
    m_scene->setModel(m_sketcher_model);
    m_ui->top_bar_wdg->setModel(m_sketcher_model);
    m_ui->side_bar_wdg->setModel(m_sketcher_model);

    connectTopBarSlots();
    connectSideBarSlots();

    setStyleSheet(schrodinger::sketcher::SKETCHER_WIDGET_STYLE);

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
    // Connect "More Actions" menu
    connect(m_ui->top_bar_wdg, &SketcherTopBar::selectAllRequested, m_scene,
            &Scene::selectAll);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::clearSelectionRequested,
            m_scene, &Scene::clearSelection);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::invertSelectionRequested,
            m_scene, &Scene::invertSelection);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::pasteRequested, m_scene,
            &Scene::onPasteRequested);

    // Clear/Import/Export
    connect(m_ui->top_bar_wdg, &SketcherTopBar::clearSketcherRequested, m_scene,
            &Scene::clearInteractiveItems);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::importTextRequested, m_scene,
            &Scene::onImportTextRequested);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::saveImageRequested, m_scene,
            &Scene::showFileSaveImageDialog);
    connect(m_ui->top_bar_wdg, &SketcherTopBar::exportToFileRequested, m_scene,
            &Scene::showFileExportDialog);
}

void SketcherWidget::connectSideBarSlots()
{
    // Connect "Select Options" widget
    connect(m_ui->side_bar_wdg, &SketcherSideBar::selectAllRequested, m_scene,
            &Scene::selectAll);
    connect(m_ui->side_bar_wdg, &SketcherSideBar::clearSelectionRequested,
            m_scene, &Scene::clearSelection);
    connect(m_ui->side_bar_wdg, &SketcherSideBar::invertSelectionRequested,
            m_scene, &Scene::invertSelection);

    // TODO: Testing components to be removed from the widget
    m_ui->font_sb->setValue(m_scene->fontSize());
    connect(m_ui->font_sb, &QDoubleSpinBox::valueChanged, m_scene,
            &Scene::setFontSize);
    m_ui->carbon_labels_combo->addItems(
        {"No carbon labels", "Terminal carbons only", "All atoms labeled"});
    connect(m_ui->carbon_labels_combo, &QComboBox::currentIndexChanged,
            [this](auto i) { m_scene->setCarbonsLabeled(CarbonLabels(i)); });
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