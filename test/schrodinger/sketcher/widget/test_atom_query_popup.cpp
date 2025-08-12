#define BOOST_TEST_MODULE Test_Sketcher
#include <QSignalSpy>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_atom_query_popup.h"
#include "schrodinger/sketcher/widget/atom_query_popup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Subclass that allows access to protected data for testing purposes.
 */
class TestAtomQueryPopup : public AtomQueryPopup
{
  public:
    TestAtomQueryPopup()
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

    TestAtomQueryPopup popup;
    QSignalSpy spy(&popup, &AtomQueryPopup::selectionChanged);

    for (auto packet : popup.getButtonPackets()) {
        packet.button->click();

        BOOST_TEST(spy.count() == 1);
        auto args = spy.takeLast();
        BOOST_TEST(args.at(0).typeId() == QMetaType::Int);
        BOOST_TEST(args.at(0).toInt() == static_cast<int>(packet.enum_int));
    }
}

bool view_synchronized_to_model(TestAtomQueryPopup& wdg)
{
    auto model = wdg.getModel();
    auto draw_tool = model->getDrawTool();
    auto atom_tool = model->getAtomTool();
    int exp_button_id = -1;
    if (draw_tool == DrawTool::ATOM && atom_tool == AtomTool::QUERY) {
        exp_button_id = model->getValueInt(ModelKey::ATOM_QUERY);
    }
    return wdg.getGroup()->checkedId() == exp_button_id;
}

/**
 * Verify proper check state assignment.
 */
BOOST_AUTO_TEST_CASE(checked_button)
{

    TestAtomQueryPopup wdg;
    auto model = wdg.getModel();

    std::vector<DrawTool> draw_tools = {DrawTool::ATOM, DrawTool::BOND,
                                        DrawTool::CHARGE, DrawTool::RING,
                                        DrawTool::ENUMERATION};
    std::vector<AtomTool> atom_tools = {AtomTool::ELEMENT, AtomTool::QUERY};
    std::vector<AtomQuery> queries = {
        AtomQuery::A, AtomQuery::AH, AtomQuery::M, AtomQuery::MH,
        AtomQuery::Q, AtomQuery::QH, AtomQuery::X, AtomQuery::XH};
    for (auto draw_tool : draw_tools) {
        model->setValue(ModelKey::DRAW_TOOL, draw_tool);
        for (auto atom_tool : atom_tools) {
            model->setValue(ModelKey::ATOM_TOOL, atom_tool);
            for (auto query : queries) {
                model->setValue(ModelKey::ATOM_QUERY, query);
                BOOST_TEST(view_synchronized_to_model(wdg));
            }
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
