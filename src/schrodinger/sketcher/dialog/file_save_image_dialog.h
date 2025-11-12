#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/widget/sketcher_view_with_wasm_outline_fix.h"
#include "schrodinger/sketcher/molviewer/constants.h"

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
class SKETCHER_API FileSaveImagePopup : public SketcherViewWithWasmOutlineFix
{
    Q_OBJECT

    friend class FileSaveImageDialog; // to access m_ui

  public:
    FileSaveImagePopup(QWidget* parent = nullptr,
                       SketcherModel* model = nullptr);
    ~FileSaveImagePopup();

    /**
     * Extract render options from the popup widget
     */
    RenderOptions getRenderOptions() const;

  signals:
    void renderOptionsChanged();

  private:
    std::unique_ptr<Ui::FileSaveImagePopup> m_ui;
    QColor m_background_color = LIGHT_BACKGROUND_COLOR;
};

/**
 * Dialog from which to save image files
 */
class SKETCHER_API FileSaveImageDialog : public FileExportDialog
{
    Q_OBJECT

  public:
    FileSaveImageDialog(SketcherModel* model, QWidget* parent = nullptr,
                        bool is_svg_enabled = true);
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
    void reject() override;
    void showEvent(QShowEvent* event) override;

    /**
     * Override this function to a no-op since this widget doesn't have
     * different states for reaction and non-reaction export.
     */
    void setIsReactionExport(bool has_reaction) override{};

    /**
     * @return serialized image from the current scene to write to disk
     */
    QByteArray getFileContent() const override;

    /**
     * @return supported extensions for the selected format in the dialog
     */
    std::vector<std::string> getValidExtensions() const override;

    /**
     * Respond to an option changing by updating the options label
     */
    void onRenderOptionsChanged();

    FileSaveImagePopup* m_options_popup = nullptr;
    QWidget* m_options_wdg = nullptr;
    std::unique_ptr<Ui::FileSaveImageWidget> m_options_wdg_ui;
    int m_image_width_at_start = 0;
    int m_image_height_at_start = 0;
    bool m_image_transparency_at_start = false;
};

} // namespace sketcher
} // namespace schrodinger
