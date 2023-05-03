#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/abstract_bond_popup.h"

namespace Ui
{
class BondQueryPopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Popup used to provide functionality choices for a "bond query" button.
 */
class SKETCHER_API BondQueryPopup : public AbstractBondPopup
{
    Q_OBJECT

  public:
    BondQueryPopup(QWidget* parent = nullptr);
    ~BondQueryPopup();

  protected:
    std::unique_ptr<Ui::BondQueryPopup> ui;
    void generateButtonPackets() override;
};

} // namespace sketcher
} // namespace schrodinger
