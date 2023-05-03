#pragma once

#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/modal_dialog.h"

namespace Ui
{
class About2DSketcher;
}

namespace schrodinger
{
namespace sketcher
{
class SKETCHER_API About2DSketcher : public ModalDialog
{
    Q_OBJECT
  public:
    About2DSketcher(QWidget* parent = nullptr);
    ~About2DSketcher();

  private:
    std::unique_ptr<Ui::About2DSketcher> m_ui;
};

} // namespace sketcher
} // namespace schrodinger
