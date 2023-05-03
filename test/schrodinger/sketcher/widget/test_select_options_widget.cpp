#define BOOST_TEST_MODULE Test_Sketcher
#include <QButtonGroup>
#include <QTest>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/Scene.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_select_options_widget.h"
#include "schrodinger/sketcher/widget/select_options_widget.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

namespace schrodinger
{
namespace sketcher
{

std::vector<SetterPacket> m_setter_packets;

/**
 * Data necessary for updating the state of the model from checkable
 * actions on this view.
 *
 * These are used when synchronizing the model to this view.
 */
std::vector<SignalPacket> m_signal_packets;

/**
 * Subclass that allows access to protected data for testing purposes.
 */
class TestSelectOptionsWidget : public SelectOptionsWidget
{
  public:
    TestSelectOptionsWidget()
    {
        setModel(new SketcherModel(this));
    }

    Ui::SelectOptionsWidget* getUI() const
    {
        return ui.get();
    }
};

/**
 * Compare the state of the `SelectOptionsWidget`'s subwidgets and model.
 *
 * @param wdg A a select options widget
 */
bool view_synchronized_to_model(TestSelectOptionsWidget& wdg)
{
    auto model = wdg.getModel();
    auto ui = wdg.getUI();
    auto draw_tool = model->getDrawTool();

    int sel_button_id = -1;
    if (draw_tool == DrawTool::SELECT) {
        sel_button_id = model->getValueInt(ModelKey::SELECTION_TOOL);
    }
    bool move_rotate_checked = draw_tool == DrawTool::MOVE_ROTATE;
    bool erase_checked = draw_tool == DrawTool::ERASE;

    return ui->select_tool_group->checkedId() == sel_button_id &&
           ui->move_rotate_btn->isChecked() == move_rotate_checked &&
           ui->erase_btn->isChecked() == erase_checked;
}

/**
 * Verify proper synchronization with model.
 */
BOOST_AUTO_TEST_CASE(synchronize)
{
    TestSelectOptionsWidget wdg;
    auto model = wdg.getModel();
    auto ui = wdg.getUI();
    sketcherScene scene;
    scene.setModel(model);
    scene.importText("C"); // require contents for the clicks to activate

    // setModel() should synchronize the view and model immediately
    model->setValue(ModelKey::DRAW_TOOL, DrawTool::SELECT);
    BOOST_TEST(view_synchronized_to_model(wdg));

    // Try interacting with subwidgets of the view. Synchronization should occur
    // automatically.
    auto group = ui->select_tool_group;
    BOOST_TEST(group->checkedId() ==
               static_cast<int>(SelectionTool::RECTANGLE));
    group->button(static_cast<int>(SelectionTool::LASSO))->click();
    BOOST_TEST(group->checkedId() == static_cast<int>(SelectionTool::LASSO));
    BOOST_TEST(view_synchronized_to_model(wdg));
    ui->move_rotate_btn->click();
    BOOST_TEST(view_synchronized_to_model(wdg));
    ui->erase_btn->click();
    BOOST_TEST(view_synchronized_to_model(wdg));
    group->button(static_cast<int>(SelectionTool::RECTANGLE))->click();
    BOOST_TEST(view_synchronized_to_model(wdg));

    // Try changing the state of the model. The view should update the state of
    // its subwidgets automatically.
    for (auto& setter_packet : wdg.getSetterPackets()) {
        int new_value = !model->getValue(setter_packet.key).toInt();
        model->setValue(setter_packet.key, new_value);
        BOOST_TEST(model->getValue(setter_packet.key).toInt() == new_value);
        BOOST_TEST(view_synchronized_to_model(wdg));
    }
}

/**
 * Verify that widgets are enabled and disabled as expected.
 */
BOOST_AUTO_TEST_CASE(updateWidgetsEnabled)
{
    TestSelectOptionsWidget wdg;
    auto model = wdg.getModel();
    auto ui = wdg.getUI();
    QGraphicsTextItem text_item;
    sketcherScene scene;
    scene.setModel(model);

    auto requires_contents_btns = {ui->select_square_btn, ui->select_lasso_btn,
                                   ui->move_rotate_btn, ui->erase_btn,
                                   ui->select_all_btn};
    auto requires_selection_btns = {ui->clear_selection_btn,
                                    ui->invert_selection_btn};

    // Empty scene
    for (auto btn : requires_contents_btns) {
        BOOST_TEST(!btn->isEnabled());
    }
    for (auto btn : requires_selection_btns) {
        BOOST_TEST(!btn->isEnabled());
    }

    scene.importText("CC");
    for (bool has_selection : {true, false}) {
        if (has_selection) {
            scene.selectAll();
        } else {
            scene.clearSelection();
        }
        for (auto btn : requires_contents_btns) {
            BOOST_TEST(btn->isEnabled());
        }
        for (auto btn : requires_selection_btns) {
            BOOST_TEST(btn->isEnabled() == has_selection);
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
