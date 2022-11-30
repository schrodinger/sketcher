#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/sketcher_view.h"

namespace Ui
{
class FileSaveImageWidget;
class FileSaveImagePopup;
} // namespace Ui

namespace schrodinger
{
namespace sketcher
{
class SketcherModel;
enum class ImageFormat;

/**
 * Popup widget containing rendering options
 */
class SKETCHER_API FileSaveImagePopup : public SketcherView
{
    Q_OBJECT

  public:
    FileSaveImagePopup(QWidget* parent = nullptr);
    ~FileSaveImagePopup();

    /**
     * Extract render options from the popup widget
     */
    RenderOptions getRenderOptions() const;

  signals:
    void renderOptionsChanged();

  private:
    /**
     * NOTE: Required to pick up stylesheet that makes the background white
     */
    void paintEvent(QPaintEvent* event) override;

    std::unique_ptr<Ui::FileSaveImagePopup> m_ui;
};

/**
 * Dialog from which to save image files
 */
class SKETCHER_API FileSaveImageDialog : public FileExportDialog
{
    Q_OBJECT

  public:
    FileSaveImageDialog(SketcherModel* model, QWidget* parent = nullptr);
    ~FileSaveImageDialog();

  signals:
    /**
     * Signal used to request serialized image from the associated scene.
     *
     * @param format requested image format
     * @param opts image rendering options
     * @return bytes export image from the current scene
     */
    QByteArray exportImageRequested(ImageFormat format,
                                    const RenderOptions& opts) const;

  private:
    /**
     * @return serialized image from the current scene to write to disk
     */
    QByteArray getFileContent() const override;

    /**
     * @return supported extensions for the selected format in the dialog
     */
    QStringList getValidExtensions() const override;

    /**
     * Respond to an option changing by updating the options label
     */
    void onRenderOptionsChanged();

    FileSaveImagePopup* m_options_popup = nullptr;
    QWidget* m_options_wdg = nullptr;
    std::unique_ptr<Ui::FileSaveImageWidget> m_options_wdg_ui;
};

} // namespace sketcher
} // namespace schrodinger
