#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>
#include "../test_common.h"
#include "schrodinger/sketcher/Scene.h"
#include "schrodinger/sketcher/widget/enumeration_tool_widget.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_enumeration_tool_widget.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Subclass that allows access to protected data for testing purposes.
 */
class TestEnumerationToolWidget : public EnumerationToolWidget
{
  public:
    TestEnumerationToolWidget()
    {
        setModel(new SketcherModel(this));
    }
    std::unique_ptr<Ui::EnumerationToolWidget>& getUI()
    {
        return ui;
    }
};

/**
 * Verify that the enumeration buttons work as expected.
 */
BOOST_AUTO_TEST_CASE(enumeration_buttons)
{

    TestEnumerationToolWidget wdg;
    auto model = wdg.getModel();
    auto& ui = wdg.getUI();

    model->setValue(ModelKey::DRAW_TOOL,
                    QVariant(static_cast<int>(DrawTool::ENUMERATION)));

    // The check state of the button should reflect the state of the model
    for (auto enum_tool :
         {EnumerationTool::NEW_RGROUP, EnumerationTool::EXISTING_RGROUP,
          EnumerationTool::ATTACHMENT_POINT, EnumerationTool::ADD_MAPPING,
          EnumerationTool::REMOVE_MAPPING, EnumerationTool::RXN_ARROW,
          EnumerationTool::RXN_PLUS}) {
        model->setValue(ModelKey::ENUMERATION_TOOL,
                        QVariant(static_cast<int>(enum_tool)));

        bool exp_checked = enum_tool == EnumerationTool::NEW_RGROUP ||
                           enum_tool == EnumerationTool::EXISTING_RGROUP;
        BOOST_TEST(ui->rgroup_btn->isChecked() == exp_checked);
        if (exp_checked) {
            int rgroup_number = ui->rgroup_btn->getEnumItem();
            int exp_rgroup_number =
                enum_tool == EnumerationTool::NEW_RGROUP
                    ? 0
                    : model->getValue(ModelKey::RGROUP_NUMBER).toInt();
            BOOST_TEST(rgroup_number == exp_rgroup_number);
        }

        exp_checked = enum_tool == EnumerationTool::ATTACHMENT_POINT;
        BOOST_TEST(ui->attachment_point_btn->isChecked() == exp_checked);
        exp_checked = enum_tool == EnumerationTool::RXN_ARROW ||
                      enum_tool == EnumerationTool::RXN_PLUS ||
                      enum_tool == EnumerationTool::ADD_MAPPING ||
                      enum_tool == EnumerationTool::REMOVE_MAPPING;
        BOOST_TEST(ui->reaction_btn->isChecked() == exp_checked);
    }

    // Assign values to the model in order to update the modular tool buttons on
    // this widget
    model->setValue(ModelKey::ENUMERATION_TOOL,
                    QVariant(static_cast<int>(EnumerationTool::RXN_ARROW)));
    model->setValue(
        ModelKey::ENUMERATION_TOOL,
        QVariant(static_cast<int>(EnumerationTool::EXISTING_RGROUP)));
    model->setValue(ModelKey::RGROUP_NUMBER, QVariant(5));

    // Clicking these buttons should modify the model
    ui->reaction_btn->click();
    BOOST_TEST(ui->reaction_btn->isChecked());
    BOOST_TEST(model->getValue(ModelKey::ENUMERATION_TOOL).toInt() ==
               static_cast<int>(EnumerationTool::RXN_ARROW));

    ui->attachment_point_btn->click();
    BOOST_TEST(ui->attachment_point_btn->isChecked());
    BOOST_TEST(model->getValue(ModelKey::ENUMERATION_TOOL).toInt() ==
               static_cast<int>(EnumerationTool::ATTACHMENT_POINT));

    ui->attachment_point_btn->click();
    BOOST_TEST(ui->attachment_point_btn->isChecked());
    BOOST_TEST(model->getValue(ModelKey::ENUMERATION_TOOL).toInt() ==
               static_cast<int>(EnumerationTool::ATTACHMENT_POINT));
    BOOST_TEST(model->getValue(ModelKey::RGROUP_NUMBER).toInt() == 5);

    // Switching to a different interaction mode or draw tool should uncheck all
    // buttons in this widget
    model->setValue(ModelKey::ENUMERATION_TOOL,
                    QVariant(static_cast<int>(EnumerationTool::NEW_RGROUP)));
    std::vector<DrawTool> draw_tools = {DrawTool::ATOM, DrawTool::BOND,
                                        DrawTool::CHARGE, DrawTool::RING,
                                        DrawTool::ENUMERATION};
    for (auto& draw_tool : draw_tools) {
        model->setValue(ModelKey::DRAW_TOOL,
                        QVariant(static_cast<int>(draw_tool)));

        bool in_enum_mode = draw_tool == DrawTool::ENUMERATION;
        for (auto& enum_tool :
             {EnumerationTool::NEW_RGROUP, EnumerationTool::EXISTING_RGROUP,
              EnumerationTool::ATTACHMENT_POINT, EnumerationTool::ADD_MAPPING,
              EnumerationTool::REMOVE_MAPPING}) {
            model->setValue(ModelKey::ENUMERATION_TOOL,
                            QVariant(static_cast<int>(enum_tool)));
            for (unsigned int rgroup_number = 0; rgroup_number < 15;
                 ++rgroup_number) {
                model->setValue(ModelKey::RGROUP_NUMBER,
                                QVariant(rgroup_number));
                bool exp_checked =
                    in_enum_mode &&
                    (enum_tool == EnumerationTool::NEW_RGROUP ||
                     (enum_tool == EnumerationTool::EXISTING_RGROUP &&
                      (rgroup_number > 0 && rgroup_number < 10)));
                BOOST_TEST(ui->rgroup_btn->isChecked() == exp_checked);
                exp_checked = in_enum_mode &&
                              enum_tool == EnumerationTool::ATTACHMENT_POINT;
                BOOST_TEST(ui->attachment_point_btn->isChecked() ==
                           exp_checked);
                exp_checked = in_enum_mode &&
                              (enum_tool == EnumerationTool::RXN_ARROW ||
                               enum_tool == EnumerationTool::RXN_PLUS ||
                               enum_tool == EnumerationTool::ADD_MAPPING ||
                               enum_tool == EnumerationTool::REMOVE_MAPPING);
                BOOST_TEST(ui->reaction_btn->isChecked() == exp_checked);
            }
        }
    }
}

/**
 * Verify that widgets are enabled and disabled as expected.
 */
BOOST_AUTO_TEST_CASE(updateWidgetsEnabled)
{
    TestEnumerationToolWidget wdg;
    auto model = wdg.getModel();
    sketcherScene scene;
    scene.setModel(model);
    scene.importText("CC");
    BOOST_TEST(wdg.isEnabled());
    scene.selectAll();
    BOOST_TEST(!wdg.isEnabled());
}

} // namespace sketcher
} // namespace schrodinger
