#pragma once
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/modular_popup.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Abstract popup class for all bond tool popups.
 */
class SKETCHER_API AbstractBondPopup : public ModularPopup
{
    Q_OBJECT

  public:
    AbstractBondPopup(QWidget* parent = nullptr);

  protected:
    int getButtonIDToCheck() override;
};

} // namespace sketcher
} // namespace schrodinger
