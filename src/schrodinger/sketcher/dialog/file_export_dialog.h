#pragma once

#include "schrodinger/sketcher/dialog/modal_dialog.h"
#include "schrodinger/sketcher/definitions.h"

namespace Ui
{
class FileExportDialog;
}

namespace schrodinger
{
namespace rdkit_extensions
{
enum class Format;
}
namespace sketcher
{

class SketcherModel;

/**
 * Dialog from which to export structure files
 */
class SKETCHER_API FileExportDialog : public ModalDialog
{
    Q_OBJECT

  public:
    FileExportDialog(SketcherModel* model, QWidget* parent = nullptr);
    ~FileExportDialog();

  signals:
    /**
     * Signal used to request serialized text from the associated scene.
     *
     * @param format requested structure format
     * @return text export text from the current scene
     */
    QString exportTextRequested(rdkit_extensions::Format format) const;

  protected:
    /**
     * @return serialized structure from the current scene to write to disk
     */
    virtual QByteArray getFileContent() const;

    /**
     * @return supported extensions for the selected format in the dialog
     */
    virtual QStringList getValidExtensions() const;

    /**
     * Exports the current scene from the filename/format specified in the UI
     */
    void exportFile();

    std::unique_ptr<Ui::FileExportDialog> m_ui;
    bool m_model_has_reaction = false;
};

} // namespace sketcher
} // namespace schrodinger
