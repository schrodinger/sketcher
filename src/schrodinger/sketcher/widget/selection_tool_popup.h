#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace Ui
{
class SelectionToolPopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Popup used to provide functionality choices for a selection tool button.
 */
class SKETCHER_API SelectionToolPopup : public ModularPopup
{
    Q_OBJECT

  public:
    SelectionToolPopup(QWidget* parent = nullptr);
    ~SelectionToolPopup();

  protected:
    std::unique_ptr<Ui::SelectionToolPopup> ui;
    void generateButtonPackets() override;
    int getButtonIDToCheck() override;
};

} // namespace sketcher
} // namespace schrodinger
