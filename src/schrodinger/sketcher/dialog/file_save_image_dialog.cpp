#include "schrodinger/sketcher/dialog/file_save_image_dialog.h"

#include <QPainter>

#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher_css_style.h"
#include "schrodinger/sketcher/ui/ui_file_export_dialog.h"
#include "schrodinger/sketcher/ui/ui_file_save_image_popup.h"
#include "schrodinger/sketcher/ui/ui_file_save_image_widget.h"

namespace schrodinger
{
namespace sketcher
{

FileSaveImagePopup::FileSaveImagePopup(QWidget* parent, SketcherModel* model) :
    SketcherViewWithWasmOutlineFix(parent)
{
    m_ui.reset(new Ui::FileSaveImagePopup());
    m_ui->setupUi(this);

    setWindowFlags(Qt::Popup);
    setStyleSheet(RENDER_OPTIONS_POPUP_STYLE);

    connect(m_ui->width_sb, &QSpinBox::valueChanged, this,
            &FileSaveImagePopup::renderOptionsChanged);
    connect(m_ui->height_sb, &QSpinBox::valueChanged, this,
            &FileSaveImagePopup::renderOptionsChanged);
#if QT_VERSION >= QT_VERSION_CHECK(6, 7, 0)
    connect(m_ui->transparent_cb, &QCheckBox::checkStateChanged, this,
            &FileSaveImagePopup::renderOptionsChanged);
#else
    connect(m_ui->transparent_cb, &QCheckBox::stateChanged, this,
            &FileSaveImagePopup::renderOptionsChanged);
#endif
    connect(model, &SketcherModel::backgroundColorChanged, this,
            [this](const QColor& color) {
                m_background_color = color;
                renderOptionsChanged();
            });
}

FileSaveImagePopup::~FileSaveImagePopup() = default;

RenderOptions FileSaveImagePopup::getRenderOptions() const
{
    RenderOptions opts;
    opts.width_height =
        QSize(m_ui->width_sb->value(), m_ui->height_sb->value());
    opts.background_color = m_ui->transparent_cb->isChecked()
                                ? Qt::transparent
                                : m_background_color;
    return opts;
}

FileSaveImageDialog::FileSaveImageDialog(SketcherModel* model, QWidget* parent,
                                         bool is_svg_enabled) :
    FileExportDialog(model, parent)
{
    setWindowTitle("Save Image");

    // Set up the options widget UI and popup
    m_options_wdg = new QWidget(this);
    m_options_wdg_ui.reset(new Ui::FileSaveImageWidget());
    m_options_wdg_ui->setupUi(m_options_wdg);
    m_options_wdg->setStyleSheet(RENDER_OPTIONS_STYLE);

    // Setup the popup widget that handles rendering options
    m_options_popup = new FileSaveImagePopup(this, model);
    m_options_wdg_ui->change_btn->setPopupDelay(0);
    m_options_wdg_ui->change_btn->setPopupWidget(m_options_popup);
    m_options_wdg_ui->change_btn->showPopupIndicator(false);
    m_options_wdg_ui->change_btn->setStyleSheet(TEXT_LINK_STYLE);

    // Add to the last row of the file/format layout, spanning all columns
    m_ui->file_format_layout->addWidget(m_options_wdg, 2, 0, 1, -1);

    m_ui->format_combo->clear(); // clear out base class initialization
    for (const auto& [format, label, exts] : get_image_export_formats()) {
        auto filter = get_filter_name(label, exts);
        if (!is_svg_enabled && format == ImageFormat::SVG) {
            continue;
        }
        m_ui->format_combo->addItem(filter, QVariant::fromValue(format));
    }

    connect(m_options_popup, &FileSaveImagePopup::renderOptionsChanged, this,
            &FileSaveImageDialog::onRenderOptionsChanged);
    onRenderOptionsChanged();
}

FileSaveImageDialog::~FileSaveImageDialog() = default;

QByteArray FileSaveImageDialog::getFileContent() const
{
    auto combo_format = m_ui->format_combo->currentData().value<ImageFormat>();
    auto opts = m_options_popup->getRenderOptions();
    return emit exportImageRequested(combo_format, opts);
}

std::vector<std::string> FileSaveImageDialog::getValidExtensions() const
{
    auto combo_format = m_ui->format_combo->currentData().value<ImageFormat>();
    return {get_image_extension(combo_format)};
}

void FileSaveImageDialog::onRenderOptionsChanged()
{
    auto opts = m_options_popup->getRenderOptions();
    QString background_color = "Solid";
    if (opts.background_color == Qt::transparent) {
        background_color = "Transparent";
    } else if (opts.background_color == LIGHT_BACKGROUND_COLOR) {
        background_color = "White";
    } else if (opts.background_color == DARK_BACKGROUND_COLOR) {
        background_color = "Dark";
    }
    QString width = QString::number(opts.width_height.width());
    QString height = QString::number(opts.width_height.height());
    m_options_wdg_ui->status_lbl->setText(background_color + " background, " +
                                          width + " x " + height + " px");
}

void FileSaveImageDialog::reject()
{
    m_options_popup->m_ui->width_sb->setValue(m_image_width_at_start);
    m_options_popup->m_ui->height_sb->setValue(m_image_height_at_start);
    m_options_popup->m_ui->transparent_cb->setChecked(
        m_image_transparency_at_start);
    FileExportDialog::reject();
}

void FileSaveImageDialog::showEvent(QShowEvent* event)
{
    ModalDialog::showEvent(event);

    // Store the current image options to restore them later
    m_image_width_at_start = m_options_popup->m_ui->width_sb->value();
    m_image_height_at_start = m_options_popup->m_ui->height_sb->value();
    m_image_transparency_at_start =
        m_options_popup->m_ui->transparent_cb->isChecked();
    FileExportDialog::showEvent(event);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/file_save_image_dialog.moc"
