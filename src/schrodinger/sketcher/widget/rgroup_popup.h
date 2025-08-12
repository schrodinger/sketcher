#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace Ui
{
class RGroupPopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Popup used to provide functionality choices for the R group button.
 */
class SKETCHER_API RGroupPopup : public ModularPopup
{
    Q_OBJECT

  public:
    RGroupPopup(QWidget* parent = nullptr);
    ~RGroupPopup();

    /**
     * @param rgroup_number An R group number
     * @return Whether this widget contains a button corresponding to the R
     * group specified by `rgroup_number`
     */
    bool isSupportedRGroup(unsigned int rgroup_number);

  protected:
    std::unique_ptr<Ui::RGroupPopup> ui;
    void generateButtonPackets() override;
    int getButtonIDToCheck() override;
};

} // namespace sketcher
} // namespace schrodinger
