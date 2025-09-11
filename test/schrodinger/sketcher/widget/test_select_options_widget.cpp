#define BOOST_TEST_MODULE Test_Sketcher
#include <QButtonGroup>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_select_options_widget.h"
#include "schrodinger/sketcher/widget/select_options_widget.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

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
    TestSelectOptionsWidget(SketcherModel* model)
    {
        setModel(model);
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

    int checked_selection_btn_id = -1;
    auto group = ui->select_tool_group;
    for (auto* btn : group->buttons()) {
        auto mod_button = dynamic_cast<ModularToolButton*>(btn);
        if (mod_button != nullptr && mod_button->isChecked()) {
            checked_selection_btn_id = mod_button->getEnumItem();
        }
    }

    return checked_selection_btn_id == sel_button_id &&
           ui->move_rotate_btn->isChecked() == move_rotate_checked &&
           ui->erase_btn->isChecked() == erase_checked;
}

/**
 * Verify proper synchronization with model.
 */
BOOST_AUTO_TEST_CASE(synchronize)
{
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    TestSelectOptionsWidget wdg(model);
    auto ui = wdg.getUI();
    // require contents for the clicks to activate
    import_mol_text(test_scene->m_mol_model, "C");

    // setModel() should synchronize the view and model immediately
    model->setValue(ModelKey::DRAW_TOOL, DrawTool::SELECT);
    BOOST_TEST(view_synchronized_to_model(wdg));

    // Try interacting with subwidgets of the view. Synchronization should occur
    // automatically.

    auto group = ui->select_tool_group;
    BOOST_TEST_REQUIRE(group != nullptr);
    BOOST_TEST_REQUIRE(!group->buttons().isEmpty());

    int checked_rectangles = 0;
    for (auto* btn : group->buttons()) {
        auto mod_button = dynamic_cast<ModularToolButton*>(btn);
        if (mod_button != nullptr &&
            mod_button->getEnumItem() ==
                static_cast<int>(SelectionTool::RECTANGLE) &&
            mod_button->isChecked()) {
            checked_rectangles++;
            break;
        }
    }
    BOOST_TEST(checked_rectangles == 1);
    ui->select_tool_btn_2->click();
    BOOST_TEST(view_synchronized_to_model(wdg));
    ui->move_rotate_btn->click();
    BOOST_TEST(view_synchronized_to_model(wdg));
    ui->erase_btn->click();
    BOOST_TEST(view_synchronized_to_model(wdg));
    ui->select_tool_btn_1->click();
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
    auto test_scene = TestScene::getScene();
    auto mol_model = test_scene->m_mol_model;
    TestSelectOptionsWidget wdg(test_scene->m_sketcher_model);
    auto ui = wdg.getUI();
    QGraphicsTextItem text_item;

    auto requires_contents_btns = {ui->move_rotate_btn, ui->erase_btn};
    auto requires_selection_btns = {ui->clear_selection_btn,
                                    ui->invert_selection_btn};

    // Empty scene
    for (auto btn : requires_contents_btns) {
        BOOST_TEST(!btn->isEnabled());
    }
    for (auto btn : requires_selection_btns) {
        BOOST_TEST(!btn->isEnabled());
    }

    import_mol_text(mol_model, "CC");
    for (bool has_selection : {true, false}) {
        if (has_selection) {
            mol_model->selectAll();
        } else {
            mol_model->clearSelection();
        }
        // select all is disabled when everything is selected
        BOOST_TEST(ui->select_all_btn->isEnabled() == !has_selection);
        for (auto btn : requires_contents_btns) {
            BOOST_TEST(btn->isEnabled());
        }
        for (auto btn : requires_selection_btns) {
            BOOST_TEST(btn->isEnabled() == has_selection);
        }
    }

    // SKETCH-2219: Ensure that clicking any of these mutually exclusive buttons
    // twice doesn't deselect it
    for (auto btn : requires_contents_btns) {
        btn->click();
        BOOST_TEST(btn->isChecked());
        btn->click();
        BOOST_TEST(btn->isChecked());
        for (auto other_btn : requires_contents_btns) {
            BOOST_TEST(other_btn->isChecked() == (btn == other_btn));
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
