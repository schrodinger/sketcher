#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/abstract_bond_popup.h"

namespace Ui
{
class StereoBondPopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Popup used to provide functionality choices for a "cis/trans bond" button.
 */
class SKETCHER_API StereoBondPopup : public AbstractBondPopup
{
    Q_OBJECT

  public:
    StereoBondPopup(QWidget* parent = nullptr);
    ~StereoBondPopup();

  protected:
    std::unique_ptr<Ui::StereoBondPopup> ui;
    void generateButtonPackets() override;
};

} // namespace sketcher
} // namespace schrodinger
