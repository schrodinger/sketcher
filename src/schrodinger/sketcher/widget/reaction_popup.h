#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace Ui
{
class ReactionPopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Popup used to provide functionality choices for the reaction button.
 */
class SKETCHER_API ReactionPopup : public ModularPopup
{
    Q_OBJECT

  public:
    ReactionPopup(QWidget* parent = nullptr);
    ~ReactionPopup();

  protected:
    std::unique_ptr<Ui::ReactionPopup> ui;
    void generateButtonPackets() override;
    int getButtonIDToCheck() override;
    void updateWidgetsEnabled() override;
};

} // namespace sketcher
} // namespace schrodinger
