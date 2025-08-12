#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_ring_tool_widget.h"
#include "schrodinger/sketcher/widget/ring_tool_widget.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

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
    TestRingToolWidget(SketcherModel* model)
    {
        setModel(model);
    }
    using RingToolWidget::ui;
};

/**
 * Verify that the ring buttons work as expected.
 */
BOOST_AUTO_TEST_CASE(ring_buttons)
{
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    auto mol_model = test_scene->m_mol_model;
    TestRingToolWidget wdg(model);
    auto& ui = wdg.ui;
    auto group = ui->group;

    import_mol_text(mol_model, "CC");
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
    mol_model->selectAll();
    BOOST_TEST(!group->button(button_id)->isEnabled());
    mol_model->clearSelection();

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
    auto test_scene = TestScene::getScene();
    auto mol_model = test_scene->m_mol_model;
    TestRingToolWidget wdg(test_scene->m_sketcher_model);
    import_mol_text(mol_model, "CC");
    BOOST_TEST(wdg.isEnabled());
    mol_model->selectAll();
    BOOST_TEST(!wdg.isEnabled());
}

} // namespace sketcher
} // namespace schrodinger
