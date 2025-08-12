#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace Ui
{
class AtomQueryPopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Popup used to provide functionality choices for an "atom query" button.
 */
class SKETCHER_API AtomQueryPopup : public ModularPopup
{
    Q_OBJECT

  public:
    AtomQueryPopup(QWidget* parent = nullptr);
    ~AtomQueryPopup();

  protected:
    std::unique_ptr<Ui::AtomQueryPopup> ui;
    void generateButtonPackets() override;
    int getButtonIDToCheck() override;
};

} // namespace sketcher
} // namespace schrodinger
