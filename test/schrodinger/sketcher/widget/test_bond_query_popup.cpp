#define BOOST_TEST_MODULE Test_Sketcher
#include <QSignalSpy>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_bond_query_popup.h"
#include "schrodinger/sketcher/widget/bond_query_popup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Subclass that allows access to protected data for testing purposes.
 */
class TestBondQueryPopup : public BondQueryPopup
{
  public:
    TestBondQueryPopup()
    {
        setModel(new SketcherModel(this));
    }
    QButtonGroup* getGroup() const
    {
        return m_group;
    }
};

/**
 * Verify proper signal emission.
 */
BOOST_AUTO_TEST_CASE(selectionChanged)
{

    BondQueryPopup popup;
    QSignalSpy spy(&popup, &BondQueryPopup::selectionChanged);

    for (auto packet : popup.getButtonPackets()) {
        packet.button->click();

        BOOST_TEST(spy.count() == 1);
        auto args = spy.takeLast();
        BOOST_TEST(args.at(0).typeId() == QMetaType::Int);
        BOOST_TEST(args.at(0).toInt() == static_cast<int>(packet.enum_int));
    }
}

bool view_synchronized_to_model(TestBondQueryPopup& wdg)
{
    auto model = wdg.getModel();
    auto group = wdg.getGroup();
    auto draw_tool = model->getDrawTool();
    int exp_button_id = -1;
    if (draw_tool == DrawTool::BOND) {
        int bond_tool_int = model->getValueInt(ModelKey::BOND_TOOL);
        for (auto button : group->buttons()) {
            if (group->id(button) == bond_tool_int) {
                exp_button_id = bond_tool_int;
                break;
            }
        }
    }
    return group->checkedId() == exp_button_id;
}

/**
 * Verify proper check state assignment.
 */
BOOST_AUTO_TEST_CASE(checked_button)
{

    TestBondQueryPopup wdg;
    auto model = wdg.getModel();

    std::vector<DrawTool> draw_tools = {DrawTool::ATOM, DrawTool::BOND,
                                        DrawTool::CHARGE, DrawTool::RING,
                                        DrawTool::ENUMERATION};
    std::vector<BondTool> bond_tools = {BondTool::ZERO,
                                        BondTool::SINGLE,
                                        BondTool::DOUBLE,
                                        BondTool::TRIPLE,
                                        BondTool::COORDINATE,
                                        BondTool::ZERO,
                                        BondTool::SINGLE_UP,
                                        BondTool::SINGLE_DOWN,
                                        BondTool::SINGLE_EITHER,
                                        BondTool::DOUBLE_EITHER,
                                        BondTool::AROMATIC,
                                        BondTool::SINGLE_OR_DOUBLE,
                                        BondTool::SINGLE_OR_AROMATIC,
                                        BondTool::DOUBLE_OR_AROMATIC,
                                        BondTool::ANY,
                                        BondTool::ATOM_CHAIN};
    for (auto draw_tool : draw_tools) {
        model->setValue(ModelKey::DRAW_TOOL, draw_tool);
        for (auto bond_tool : bond_tools) {
            model->setValue(ModelKey::BOND_TOOL, bond_tool);
            BOOST_TEST(view_synchronized_to_model(wdg));
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
