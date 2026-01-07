#pragma once

#include <memory>
#include <optional>
#include <unordered_map>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/modal_dialog.h"

class QShowEvent;
class QString;

namespace Ui
{
class PasteInTextDialog;
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
 * Dialog with text edit field where user can type/paste structure text.
 */
class SKETCHER_API PasteInTextDialog : public ModalDialog
{
    Q_OBJECT
  public:
    PasteInTextDialog(SketcherModel* model, QWidget* parent = nullptr);
    ~PasteInTextDialog();

    /**
     * Override this method to emit `textAccepted()`.
     */
    void accept() override;

  protected:
    SketcherModel* m_sketcher_model = nullptr;
    std::unique_ptr<Ui::PasteInTextDialog> ui;
    int m_autodetect_index = 0;
    QString m_autodetect_base_name;
    /**
     * The autodetected format. Will be set to std::nullopt whenever the text
     * edit is empty or when we can't determine the format.
     */
    std::optional<rdkit_extensions::Format> m_autodetected_format =
        std::nullopt;
    /**
     * A map of all formats currently loaded into the Format combo box to their
     * name
     */
    std::unordered_map<rdkit_extensions::Format, QString> m_formats;

  protected slots:
    /**
     * Update this dialog in response to changes in the model.
     */
    void onModelValuesChanged();

    /**
     * Update the autodetection detected format
     */
    void autodetectCurrentFormat();

  signals:
    void textAccepted(const std::string& text,
                      const rdkit_extensions::Format format);
};

} // namespace sketcher
} // namespace schrodinger
