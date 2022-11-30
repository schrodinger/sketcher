#include "schrodinger/sketcher/sketcher_widget.h"

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
    ui.reset(new Ui::SketcherWidgetForm());
    ui->setupUi(this);

    m_scene = new Scene(this);
    ui->view->setScene(m_scene);

    m_sketcher_model = new SketcherModel(this);
    m_scene->setModel(m_sketcher_model);
    ui->top_bar_wdg->setModel(m_sketcher_model);
    ui->side_bar_wdg->setModel(m_sketcher_model);

    connectTopBarSlots();
    connectSideBarSlots();

    setStyleSheet(schrodinger::sketcher::SKETCHER_WIDGET_STYLE);
}

SketcherWidget::~SketcherWidget() = default;

void SketcherWidget::connectTopBarSlots()
{
    // Connect "More Actions" menu
    connect(ui->top_bar_wdg, &SketcherTopBar::pasteRequested, m_scene,
            &Scene::onPasteRequested);

    // Clear/Import/Export
    connect(ui->top_bar_wdg, &SketcherTopBar::clearSketcherRequested, m_scene,
            &Scene::clear);
    connect(ui->top_bar_wdg, &SketcherTopBar::importTextRequested, m_scene,
            &Scene::onImportTextRequested);
    connect(ui->top_bar_wdg, &SketcherTopBar::saveImageRequested, m_scene,
            &Scene::showFileSaveImageDialog);
    connect(ui->top_bar_wdg, &SketcherTopBar::exportToFileRequested, m_scene,
            &Scene::showFileExportDialog);
}

void SketcherWidget::connectSideBarSlots()
{
    // TODO: Testing components to be removed from the widget
    ui->font_sb->setValue(m_scene->fontSize());
    connect(ui->font_sb, &QDoubleSpinBox::valueChanged, m_scene,
            &Scene::setFontSize);
    ui->carbon_labels_combo->addItems(
        {"No carbon labels", "Terminal carbons only", "All atoms labeled"});
    connect(ui->carbon_labels_combo, &QComboBox::currentIndexChanged,
            [this](auto i) { m_scene->setCarbonsLabeled(CarbonLabels(i)); });
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/sketcher_widget.moc"