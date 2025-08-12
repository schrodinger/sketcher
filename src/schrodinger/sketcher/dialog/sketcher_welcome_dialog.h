#pragma once

#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/modal_dialog.h"

namespace Ui
{
class SketcherWelcomeDialog;
}

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API SketcherWelcomeDialog : public ModalDialog
{
    Q_OBJECT

  public:
    SketcherWelcomeDialog(QWidget* parent = nullptr);
    ~SketcherWelcomeDialog();

    /**
     * @param visible whether to show the "Do not show again" checkbox
     */
    void setDoNotShowCheckboxVisible(bool visible);

  signals:
    /**
     *  Notify client when 'Do not show again' checkbox state is changed.
     *  @param state is the current state of checkbox.
     */
    void doNotShowCheckboxStateChanged(int state);

  private:
    std::unique_ptr<Ui::SketcherWelcomeDialog> m_ui;
};

} // namespace sketcher
} // namespace schrodinger
