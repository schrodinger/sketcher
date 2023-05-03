#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/Scene.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_ring_tool_widget.h"
#include "schrodinger/sketcher/widget/ring_tool_widget.h"

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Subclass that allows access to protected data for testing purposes.
 */
class TestRingToolWidget : public RingToolWidget
{
  public:
    TestRingToolWidget()
    {
        setModel(new SketcherModel(this));
    }
    std::unique_ptr<Ui::RingToolWidget>& getUI()
    {
        return ui;
    }
};

/**
 * Verify that the ring buttons work as expected.
 */
BOOST_AUTO_TEST_CASE(ring_buttons)
{

    TestRingToolWidget wdg;
    auto model = wdg.getModel();
    auto& ui = wdg.getUI();
    auto group = ui->group;
    sketcherScene scene;
    scene.setModel(model);
    scene.importText("CC");
    model->setValue(ModelKey::DRAW_TOOL, DrawTool::RING);

    for (auto button : group->buttons()) {
        auto button_id = group->id(button);
        button->click();
        BOOST_TEST(button->isChecked());
        BOOST_TEST(group->checkedId() == button_id);
        BOOST_TEST(model->getValueInt(ModelKey::RING_TOOL) == button_id);
    }

    // Switching to a different interaction mode should uncheck all ring buttons
    int button_id = model->getValueInt(ModelKey::RING_TOOL);
    BOOST_TEST(group->button(button_id)->isEnabled());
    scene.selectAll();
    BOOST_TEST(!group->button(button_id)->isEnabled());
    scene.clearSelection();

    // Finally, try changing the draw tool. If not set to ring, no buttons
    // should be selected
    for (auto& draw_tool : {DrawTool::ATOM, DrawTool::BOND, DrawTool::CHARGE,
                            DrawTool::RING, DrawTool::ENUMERATION}) {
        model->setValue(ModelKey::DRAW_TOOL, draw_tool);
        int exp_id = draw_tool != DrawTool::RING
                         ? -1
                         : model->getValueInt(ModelKey::RING_TOOL);
        BOOST_TEST(group->checkedId() == exp_id);
    }
}

/**
 * Verify that widgets are enabled and disabled as expected.
 */
BOOST_AUTO_TEST_CASE(updateWidgetsEnabled)
{
    TestRingToolWidget wdg;
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
