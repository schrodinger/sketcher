#pragma once
#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/abstract_draw_tool_widget.h"

namespace Ui
{
class DrawToolsWidget;
}

namespace schrodinger
{
namespace sketcher
{

enum class BondTool;
class BondOrderPopup;
class BondQueryPopup;
class StereoBondPopup;

/**
 * Interface for user to select various tools, including atoms and bonds.
 */
class SKETCHER_API DrawToolsWidget : public AbstractDrawToolWidget
{
    Q_OBJECT

  public:
    DrawToolsWidget(QWidget* parent = nullptr);
    ~DrawToolsWidget();
    void setModel(SketcherModel* model) override;

  protected:
    std::unique_ptr<Ui::DrawToolsWidget> ui;
    BondOrderPopup* m_bond_order_wdg = nullptr;
    BondQueryPopup* m_bond_query_wdg = nullptr;
    StereoBondPopup* m_stereo_bond1_wdg = nullptr;
    StereoBondPopup* m_stereo_bond2_wdg = nullptr;

    void connectLocalSlots() override;
    void updateWidgetsEnabled() override;
    void updateCheckedButton() override;
    std::unordered_set<QAbstractButton*> getCheckableButtons() override;
    void setCheckedBondButton(int bond_item);
    void setChargeTool(int charge_tool_int);
    void setDrawTool(int draw_tool_int);

    /**
     * Get the bond button that corresponds to the specified button tool.
     *
     * If a modular button can be equipped with the specified tool but is not
     * currently, it will be assigned the tool and returned only if no other
     * button is currently equipped with the specified tool. If multiple
     * modular buttons are equipped with the same tool, one of them will be
     * chosen arbitrarily.
     *
     * @param tool A bond tool value
     * @return A button corresponding to the specified value
     */
    QAbstractButton* getBondButton(BondTool tool);

  protected slots:
    void onBondButtonClicked(QAbstractButton* button);
    void onChargeButtonClicked(int button_id);
    void onExplicitHButtonClicked();
};

} // namespace sketcher
} // namespace schrodinger
