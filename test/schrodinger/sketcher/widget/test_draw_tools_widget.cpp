#define BOOST_TEST_MODULE Test_Sketcher
#include <algorithm>

#include <QSignalSpy>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_draw_tools_widget.h"
#include "schrodinger/sketcher/widget/bond_order_popup.h"
#include "schrodinger/sketcher/widget/bond_query_popup.h"
#include "schrodinger/sketcher/widget/draw_tools_widget.h"
#include "schrodinger/sketcher/widget/stereo_bond_popup.h"

Q_DECLARE_METATYPE(schrodinger::sketcher::ModelKey);
Q_DECLARE_METATYPE(std::unordered_set<schrodinger::sketcher::ModelKey>);
Q_DECLARE_METATYPE(std::unordered_set<QGraphicsItem*>);

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Subclass that allows access to protected data for testing purposes.
 */
class TestDrawToolsWidget : public DrawToolsWidget
{
  public:
    TestDrawToolsWidget(SketcherModel* model)
    {
        setModel(model);
        m_bond_order_wdg->setAttribute(Qt::WA_DontShowOnScreen, true);
        m_bond_query_wdg->setAttribute(Qt::WA_DontShowOnScreen, true);
        m_stereo_bond1_wdg->setAttribute(Qt::WA_DontShowOnScreen, true);
        m_stereo_bond2_wdg->setAttribute(Qt::WA_DontShowOnScreen, true);
    }

    using DrawToolsWidget::getBondButton;
    using DrawToolsWidget::getCheckableButtons;
    using DrawToolsWidget::ui;

  public slots:
    using DrawToolsWidget::onBondButtonClicked;
};

/**
 * Verify that the bond buttons work as expected.
 */
BOOST_AUTO_TEST_CASE(bond_buttons)
{
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    TestDrawToolsWidget wdg(model);
    auto& ui = wdg.ui;
    import_mol_text(test_scene->m_mol_model, "CC");
    qRegisterMetaType<ModelKey>("ModelKey");
    qRegisterMetaType<std::unordered_set<ModelKey>>(
        "std::unordered_set<ModelKey>");
    QSignalSpy pinged_spy(model, &SketcherModel::valuePinged);
    QSignalSpy changed_spy(model, &SketcherModel::valuesChanged);

    // First, test the non-modal buttons (because they are easy)
    pinged_spy.clear();
    ui->single_bond_btn->click();
    int exp_value_int = static_cast<int>(BondTool::SINGLE);
    BOOST_TEST(model->getValueInt(ModelKey::BOND_TOOL) == exp_value_int);
    BOOST_TEST(pinged_spy.count() == 2);

    pinged_spy.clear();
    ui->atom_chain_btn->click();
    exp_value_int = static_cast<int>(BondTool::ATOM_CHAIN);
    BOOST_TEST(model->getValueInt(ModelKey::BOND_TOOL) ==
               static_cast<int>(BondTool::ATOM_CHAIN));
    BOOST_TEST(pinged_spy.count() == 2);

    std::vector<BondTool> bond_query_tools = {
        BondTool::ANY,
        BondTool::SINGLE_OR_DOUBLE,
        BondTool::SINGLE_OR_AROMATIC,
        BondTool::DOUBLE_OR_AROMATIC,
    };

    std::vector<BondTool> bond_order_tools = {
        BondTool::DOUBLE,
        BondTool::TRIPLE,
        BondTool::COORDINATE,
        BondTool::ZERO,
    };

    std::vector<BondTool> bond_stereo_tools = {
        BondTool::SINGLE_UP,
        BondTool::SINGLE_DOWN,
        BondTool::SINGLE_EITHER,
        BondTool::DOUBLE_EITHER,
    };

    std::vector<std::pair<ModularToolButton*, std::vector<BondTool>>>
        button_tools_pairs = {
            std::make_pair(ui->bond_query_btn, bond_query_tools),
            std::make_pair(ui->bond_order_btn, bond_order_tools),
            std::make_pair(ui->stereo_bond1_btn, bond_stereo_tools),
            std::make_pair(ui->stereo_bond2_btn, bond_stereo_tools),
        };

    for (auto& pair : button_tools_pairs) {
        auto button = pair.first;
        for (auto& tool : pair.second) {
            pinged_spy.clear();
            auto int_value = static_cast<int>(tool);
            button->setEnumItem(int_value);
            BOOST_TEST((button->toolTip()).contains("press & hold") == true);
            BOOST_TEST(model->getValueInt(ModelKey::BOND_TOOL) != int_value);
            wdg.onBondButtonClicked(button);
            BOOST_TEST(model->getValueInt(ModelKey::BOND_TOOL) == int_value);
            BOOST_TEST(pinged_spy.count() == 2);
        }
    }

    // Try changing the draw tool -- if not set to BOND, this should uncheck all
    // of the buttons
    std::vector<DrawTool> draw_tools = {DrawTool::ATOM, DrawTool::BOND,
                                        DrawTool::CHARGE, DrawTool::RING,
                                        DrawTool::ENUMERATION};
    for (auto& draw_tool : draw_tools) {
        model->setValue(ModelKey::DRAW_TOOL, draw_tool);
        int bond_tool_int = -1;
        if (draw_tool == DrawTool::BOND) {
            bond_tool_int = model->getValueInt(ModelKey::BOND_TOOL);
        }

        int button_value = -1;
        QAbstractButton* button = ui->bond_group->checkedButton();
        auto button_cast = dynamic_cast<ModularToolButton*>(button);
        if (button_cast == nullptr) {
            button_value = ui->bond_group->id(button);
        } else {
            // This is a modular button, so it has a single button ID but
            // can represent many enum values.
            button_value = static_cast<int>(button_cast->getEnumItem());
        }
        BOOST_TEST(button_value == bond_tool_int);
    }

    // Change the interaction mode; buttons should become unchecked and
    // uncheckable outside of DRAW mode

    // First, check a bond button so we have something to test
    ui->single_bond_btn->click();

    auto is_checked = [](QAbstractButton* btn) { return btn->isChecked(); };
    auto is_checkable = [](QAbstractButton* btn) { return btn->isCheckable(); };
    auto btns = wdg.getCheckableButtons();
    if (!model->hasActiveSelection()) {
        BOOST_TEST(std::any_of(btns.begin(), btns.end(), is_checked));
        BOOST_TEST(std::all_of(btns.begin(), btns.end(), is_checkable));
    } else {
        BOOST_TEST(std::none_of(btns.begin(), btns.end(), is_checked));
        BOOST_TEST(std::none_of(btns.begin(), btns.end(), is_checkable));
    }

    // If an active selection, clicking buttons should not change the model
    test_scene->m_mol_model->selectAll();
    pinged_spy.clear();
    changed_spy.clear();

    int exp_model_value_int = model->getValueInt(ModelKey::BOND_TOOL);
    int exp_key_int = static_cast<int>(ModelKey::BOND_TOOL);
    for (auto& pair : button_tools_pairs) {
        auto button = pair.first;
        for (auto& tool : pair.second) {
            auto exp_value_int = static_cast<int>(tool);
            button->setEnumItem(exp_value_int);
            wdg.onBondButtonClicked(button);
            BOOST_TEST(pinged_spy.count() == 1);
            BOOST_TEST(changed_spy.count() == 0);
            auto args = pinged_spy.takeLast();
            auto key_int = static_cast<int>(args.at(0).value<ModelKey>());
            auto value_int = args.at(1).value<QVariant>().toInt();
            BOOST_TEST(key_int == exp_key_int);
            BOOST_TEST(value_int == exp_value_int);
            BOOST_TEST(model->getValueInt(ModelKey::BOND_TOOL) ==
                       exp_model_value_int);
        }
    }
}

/**
 * Verify that the other (charge & explicit H) buttons work as expected.
 */
BOOST_AUTO_TEST_CASE(other_buttons)
{
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    TestDrawToolsWidget wdg(model);
    auto& ui = wdg.ui;
    auto group = ui->charge_group;

    auto buttons = {ui->increase_charge_btn, ui->decrease_charge_btn,
                    ui->explicit_h_btn};
    for (auto button : buttons) {
        auto button_id = group->id(button);
        button->click();
        BOOST_TEST(button->isChecked());
        BOOST_TEST(group->checkedId() == button_id);
        BOOST_TEST(ui->explicit_h_btn->isChecked() ==
                   (button == ui->explicit_h_btn));
        if (button == ui->explicit_h_btn) {
            BOOST_TEST(model->getDrawTool() == DrawTool::EXPLICIT_H);
        } else {
            BOOST_TEST(model->getValueInt(ModelKey::CHARGE_TOOL) == button_id);
        }
    }

    // Finally, try changing the draw tool/fragment type. If not set to
    // CHARGE or EXPLICIT_H, nothing should be checked
    std::array<DrawTool, 6> draw_tools = {
        DrawTool::ATOM, DrawTool::BOND,        DrawTool::CHARGE,
        DrawTool::RING, DrawTool::ENUMERATION, DrawTool::EXPLICIT_H};
    for (auto& draw_tool : draw_tools) {
        model->setValue(ModelKey::DRAW_TOOL, draw_tool);
        int exp_id = -1;
        if (draw_tool == DrawTool::CHARGE) {
            exp_id = model->getValueInt(ModelKey::CHARGE_TOOL);
        }
        BOOST_TEST(group->checkedId() == exp_id);
        BOOST_TEST(ui->explicit_h_btn->isChecked() ==
                   (draw_tool == DrawTool::EXPLICIT_H));
    }

    // SKETCH-2219: Ensure that clicking any of these mutually exclusive buttons
    // twice doesn't deselect it
    for (auto btn : buttons) {
        btn->click();
        BOOST_TEST(btn->isChecked());
        btn->click();
        BOOST_TEST(btn->isChecked());
        for (auto other_btn : buttons) {
            BOOST_TEST(other_btn->isChecked() == (btn == other_btn));
        }
    }
}

/**
 * Verify that widgets are enabled and disabled as expected.
 */
BOOST_AUTO_TEST_CASE(updateWidgetsEnabled)
{
    auto test_scene = TestScene::getScene();
    auto mol_model = test_scene->m_mol_model;
    TestDrawToolsWidget wdg(test_scene->m_sketcher_model);
    import_mol_text(mol_model, "C** |$;;_R1$|");

    auto atom = mol_model->getMol()->getAtomWithIdx(0);
    auto query = mol_model->getMol()->getAtomWithIdx(1);
    auto rgroup = mol_model->getMol()->getAtomWithIdx(2);
    auto bond = mol_model->getMol()->getBondWithIdx(0);

    auto& ui = wdg.ui;
    std::array<QWidget*, 3> other_widgets = {
        ui->increase_charge_btn, ui->decrease_charge_btn, ui->explicit_h_btn};
    std::array<QWidget*, 5> bond_widgets = {
        ui->single_bond_btn, ui->bond_order_btn, ui->bond_query_btn,
        ui->stereo_bond1_btn, ui->stereo_bond2_btn};

    auto test_enabled = [](auto widgets, auto enabled) {
        BOOST_TEST(std::all_of(
            widgets.cbegin(), widgets.cend(),
            [enabled](auto widget) { return widget->isEnabled() == enabled; }));
    };

    // No selection
    mol_model->clearSelection();
    test_enabled(other_widgets, true);
    test_enabled(bond_widgets, true);

    // Real atoms
    mol_model->select({atom}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    test_enabled(other_widgets, true);
    test_enabled(bond_widgets, false);

    // Non-atoms
    mol_model->select({query, rgroup}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    test_enabled(other_widgets, true);
    test_enabled(bond_widgets, false);

    // All atom objects
    mol_model->select({atom, query, rgroup}, {}, {}, {}, {},
                      SelectMode::SELECT_ONLY);
    test_enabled(other_widgets, true);
    test_enabled(bond_widgets, false);

    // Bond
    mol_model->select({}, {bond}, {}, {}, {}, SelectMode::SELECT_ONLY);
    test_enabled(other_widgets, false);
    test_enabled(bond_widgets, true);

    // Everything
    mol_model->selectAll();
    test_enabled(other_widgets, true);
    test_enabled(bond_widgets, true);
}

/**
 * Verify that the widget knows the correct bond button for each tool, and that
 * if the model is changed to indicate a tool not currently equipped to a
 * modular tool button, that button changes to show the tool.
 */
BOOST_AUTO_TEST_CASE(getBondButton)
{
    auto test_scene = TestScene::getScene();
    TestDrawToolsWidget wdg(test_scene->m_sketcher_model);
    auto& ui = wdg.ui;

    std::unordered_map<BondTool, QAbstractButton*> bondtool_btn_map = {
        {BondTool::SINGLE, ui->single_bond_btn},
        {BondTool::DOUBLE, ui->bond_order_btn},
        {BondTool::TRIPLE, ui->bond_order_btn},
        {BondTool::COORDINATE, ui->bond_order_btn},
        {BondTool::ZERO, ui->bond_order_btn},
        {BondTool::SINGLE_UP, ui->stereo_bond1_btn},
        {BondTool::SINGLE_DOWN, ui->stereo_bond2_btn},
        {BondTool::SINGLE_EITHER, ui->stereo_bond1_btn},
        {BondTool::DOUBLE_EITHER, ui->stereo_bond1_btn},
        {BondTool::AROMATIC, ui->bond_query_btn},
        {BondTool::SINGLE_OR_DOUBLE, ui->bond_query_btn},
        {BondTool::SINGLE_OR_AROMATIC, ui->bond_query_btn},
        {BondTool::DOUBLE_OR_AROMATIC, ui->bond_query_btn},
        {BondTool::ANY, ui->bond_query_btn},
        {BondTool::ATOM_CHAIN, ui->atom_chain_btn},
    };

    for (auto pair : bondtool_btn_map) {
        auto tool = pair.first;
        auto exp_btn = pair.second;
        auto btn = wdg.getBondButton(tool);
        BOOST_TEST(btn == exp_btn);
        auto mod_btn = dynamic_cast<ModularToolButton*>(btn);
        if (mod_btn == nullptr) {
            BOOST_TEST(ui->bond_group->id(btn) == static_cast<int>(tool));
        } else {
            BOOST_TEST(mod_btn->getEnumItem() == static_cast<int>(tool));
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
