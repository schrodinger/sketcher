#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/abstract_bond_popup.h"

namespace Ui
{
class BondOrderPopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Popup used to provide functionality choices for a "bond order" button.
 */
class SKETCHER_API BondOrderPopup : public AbstractBondPopup
{
    Q_OBJECT

  public:
    BondOrderPopup(QWidget* parent = nullptr);
    ~BondOrderPopup();

  protected:
    std::unique_ptr<Ui::BondOrderPopup> ui;
    void generateButtonPackets() override;
};

} // namespace sketcher
} // namespace schrodinger
