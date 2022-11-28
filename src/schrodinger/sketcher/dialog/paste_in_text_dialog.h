#pragma once

#include <memory>

#include "schrodinger/sketcher/dialog/modal_dialog.h"
#include "schrodinger/sketcher/definitions.h"

class QString;

namespace Ui
{
class PasteInTextDialog;
}

namespace schrodinger
{
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

  protected slots:
    /**
     * Update this dialog in response to changes in the model.
     */
    void onModelValuesChanged();

  signals:
    void textAccepted(const QString& text);
};

} // namespace sketcher
} // namespace schrodinger
