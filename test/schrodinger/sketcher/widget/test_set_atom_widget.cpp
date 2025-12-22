#define BOOST_TEST_MODULE Test_Sketcher
#include <QSignalSpy>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/ui/ui_set_atom_widget.h"
#include "schrodinger/sketcher/widget/atom_query_popup.h"
#include "schrodinger/sketcher/widget/set_atom_widget.h"

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
class TestSetAtomWidget : public SetAtomWidget
{
  public:
    TestSetAtomWidget(SketcherModel* model)
    {
        setModel(model);
        m_atom_query_wdg->setAttribute(Qt::WA_DontShowOnScreen, true);
    }
    using SetAtomWidget::getElementForButton;
    using SetAtomWidget::m_button_element_bimap;
    using SetAtomWidget::ui;

  public slots:
    using SetAtomWidget::onAtomButtonClicked;
};

/**
 * Compare the state of the `SelectOptionsWidget`'s subwidgets and model.
 *
 * @param wdg A set atom widget
 */
bool view_is_synchronized_to_model(TestSetAtomWidget& wdg)
{
    auto model = wdg.getModel();
    auto& ui = wdg.ui;
    QAbstractButton* exp_button = nullptr;
    auto draw_tool = model->getDrawTool();
    if (draw_tool == DrawTool::ATOM) {
        // The "atom" draw tool is selected, and we must dig deeper
        auto atom_tool = model->getAtomTool();
        if (atom_tool == AtomTool::ELEMENT) {
            auto bimap = wdg.m_button_element_bimap;
            auto element = model->getElement();
            if (bimap.right.count(element) == 1) {
                exp_button = bimap.right.at(element);
            } else {
                exp_button = ui->last_picked_element_btn;
                if (ui->last_picked_element_btn->getElement() != element) {
                    return false;
                }
            }
        } else {
            auto query = model->getAtomQuery();
            if (ui->atom_query_btn->getEnumItem() != static_cast<int>(query)) {
                return false;
            }
            exp_button = ui->atom_query_btn;
        }
    }
    return ui->atom_group->checkedButton() == exp_button;
}

/**
 * Verify that the view and model behave appropriately to button checking.
 */
BOOST_AUTO_TEST_CASE(button_check_state)
{
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    TestSetAtomWidget wdg(model);

    // setModel() should synchronize the view and model immediately
    BOOST_TEST(view_is_synchronized_to_model(wdg));

    // Try interacting with subwidgets of the view. Synchronization should occur
    // automatically.
    auto group = wdg.ui->atom_group;
    for (auto button : group->buttons()) {
        button->click();
        BOOST_TEST(view_is_synchronized_to_model(wdg));
        BOOST_TEST(button->isChecked());
    }

    // Try changing the state of the model. The view should update the state of
    // its subwidgets automatically.
    std::vector<AtomTool> atom_tools = {
        AtomTool::ELEMENT,
        AtomTool::QUERY,
    };
    std::vector<Element> elements = {
        Element::C, Element::H, Element::N,  Element::O,  Element::P,
        Element::S, Element::F, Element::RH, Element::SI,
    };
    std::unordered_set<AtomQuery> queries = {
        AtomQuery::A, AtomQuery::AH, AtomQuery::M, AtomQuery::MH,
        AtomQuery::Q, AtomQuery::QH, AtomQuery::X, AtomQuery::XH};
    for (auto& atom_tool : atom_tools) {
        model->setValue(ModelKey::ATOM_TOOL, atom_tool);
        BOOST_TEST(model->getAtomTool() == atom_tool);
        for (auto& element : elements) {
            model->setValue(ModelKey::ELEMENT, element);
            BOOST_TEST(model->getElement() == element);
            BOOST_TEST(view_is_synchronized_to_model(wdg));
        }
        for (auto& query : queries) {
            model->setValue(ModelKey::ATOM_QUERY, query);
            BOOST_TEST(model->getAtomQuery() == query);
            BOOST_TEST(view_is_synchronized_to_model(wdg));
        }
    }

    // Finally, try changing the draw tool -- if not set to ATOM, this should
    // uncheck all of the buttons
    std::vector<DrawTool> draw_tools = {DrawTool::ATOM, DrawTool::BOND,
                                        DrawTool::CHARGE, DrawTool::RING,
                                        DrawTool::ENUMERATION};
    for (auto& draw_tool : draw_tools) {
        model->setValue(ModelKey::DRAW_TOOL, draw_tool);
        BOOST_TEST(view_is_synchronized_to_model(wdg));
    }
}

/**
 * Verify that the modular "atom query" button works as expected.
 */
BOOST_AUTO_TEST_CASE(modular_button)
{
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    TestSetAtomWidget wdg(model);

    std::vector<AtomQuery> atom_queries = {
        AtomQuery::AH, AtomQuery::A,  AtomQuery::Q, AtomQuery::QH,
        AtomQuery::M,  AtomQuery::MH, AtomQuery::X, AtomQuery::XH,
    };

    auto button = wdg.ui->atom_query_btn;

    for (auto& atom_query : atom_queries) {
        auto int_value = static_cast<int>(atom_query);
        button->setEnumItem(int_value);
        BOOST_TEST((button->toolTip()).contains("press & hold") == true);
        BOOST_TEST(model->getAtomQuery() != atom_query);
        wdg.onAtomButtonClicked(button);
        BOOST_TEST(model->getAtomQuery() == atom_query);
    }
}

/**
 * Verify that the widget is appropriately enabled.
 */
BOOST_AUTO_TEST_CASE(enabled)
{
    auto test_scene = TestScene::getScene();
    auto mol_model = test_scene->m_mol_model;
    TestSetAtomWidget wdg(test_scene->m_sketcher_model);
    import_mol_text(mol_model, "CC");

    auto atom = mol_model->getMol()->getAtomWithIdx(0);
    auto bond = mol_model->getMol()->getBondWithIdx(0);

    // Empty selection
    mol_model->select({}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(wdg.isEnabled());
    // Atom-only selection
    mol_model->select({atom}, {}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(wdg.isEnabled());
    // Bond-only selection
    mol_model->select({}, {bond}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(!wdg.isEnabled());
    // Both atoms and bonds selected
    mol_model->select({atom}, {bond}, {}, {}, {}, SelectMode::SELECT_ONLY);
    BOOST_TEST(wdg.isEnabled());
}

/**
 * Verify that the model emits appropriate signals when the user interacts with
 * the widget.
 */
BOOST_AUTO_TEST_CASE(model_signals)
{
    auto test_scene = TestScene::getScene();
    auto model = test_scene->m_sketcher_model;
    TestSetAtomWidget wdg(model);
    import_mol_text(test_scene->m_mol_model, "CC");

    qRegisterMetaType<ModelKey>("ModelKey");
    qRegisterMetaType<std::unordered_set<ModelKey>>(
        "std::unordered_set<ModelKey>");

    QSignalSpy pinged_spy(model, &SketcherModel::valuePinged);
    QSignalSpy changed_spy(model, &SketcherModel::valuesChanged);

    // If the user clicks a button, we expect the "changed" signals to be
    // grouped together and the "set" signals to be separate.
    wdg.ui->f_btn->click();
    BOOST_TEST(changed_spy.count() == 1);
    BOOST_TEST(pinged_spy.count() == 3);

    std::unordered_set<ModelKey> pinged_keys;
    std::unordered_set<ModelKey> changed_keys;
    for (auto iter = pinged_spy.begin(); iter != pinged_spy.end(); ++iter) {
        auto pinged_args = *iter;
        pinged_keys.insert(pinged_args.at(0).value<ModelKey>());
    }

    auto changed_args = changed_spy.takeLast();
    changed_keys = changed_args.at(0).value<std::unordered_set<ModelKey>>();

    std::function<bool(ModelKey)> is_pinged_key = [pinged_keys](ModelKey key) {
        return pinged_keys.count(key) == 1;
    };

    // "Pinged" signals are the same every time
    std::unordered_set<ModelKey> exp_keys = {
        ModelKey::DRAW_TOOL, ModelKey::ATOM_TOOL, ModelKey::ELEMENT};
    BOOST_TEST(pinged_keys == exp_keys);

    // "Changed" signals depend on the current state of the model, but will
    // always affect a subset of the keys in the "set" signals
    BOOST_TEST(
        std::all_of(changed_keys.begin(), changed_keys.end(), is_pinged_key));

    // Now try with an atom query button
    pinged_spy.clear();
    wdg.ui->atom_query_btn->click();
    BOOST_TEST(changed_spy.count() == 1);
    BOOST_TEST(pinged_spy.count() == 3);

    pinged_keys.clear();
    for (auto iter = pinged_spy.begin(); iter != pinged_spy.end(); ++iter) {
        auto pinged_args = *iter;
        pinged_keys.insert(pinged_args.at(0).value<ModelKey>());
    }

    changed_args = changed_spy.takeLast();
    changed_keys = changed_args.at(0).value<std::unordered_set<ModelKey>>();

    // "Pinged" signals are the same every time
    exp_keys = {ModelKey::DRAW_TOOL, ModelKey::ATOM_TOOL, ModelKey::ATOM_QUERY};
    BOOST_TEST(pinged_keys == exp_keys);

    // "Changed" signals depend on the current state of the model, but will
    // always affect a subset of the keys in the "set" signals
    BOOST_TEST(
        std::all_of(changed_keys.begin(), changed_keys.end(), is_pinged_key));

    // Test behavior with active selection so that buttons are enabled
    auto& ui = wdg.ui;
    test_scene->m_mol_model->selectAll();
    pinged_spy.clear();
    changed_spy.clear();
    auto model_element = model->getElement();
    for (auto button : ui->atom_group->buttons()) {
        button->click();
        int exp_value_int;
        int exp_key_int;
        if (button == ui->atom_query_btn) {
            auto mod_button = dynamic_cast<ModularToolButton*>(button);
            exp_key_int = static_cast<int>(ModelKey::ATOM_QUERY);
            exp_value_int = static_cast<int>(mod_button->getEnumItem());
        } else {
            exp_key_int = static_cast<int>(ModelKey::ELEMENT);
            exp_value_int = static_cast<int>(wdg.getElementForButton(button));
        }
        BOOST_TEST(pinged_spy.count() == 1);
        BOOST_TEST(changed_spy.count() == 0);
        auto args = pinged_spy.takeLast();
        auto key_int = static_cast<int>(args.at(0).value<ModelKey>());
        auto value = args.at(1).value<QVariant>();
        BOOST_TEST(key_int == exp_key_int);
        BOOST_TEST(value.toInt() == exp_value_int);
        BOOST_TEST(model->getElement() == model_element);
    }
}

/**
 * Verify that the last-picked element button only shows elements not available
 * from the other buttons in the palette.
 */
BOOST_AUTO_TEST_CASE(last_picked_element_button)
{
    SketcherModel model;
    TestSetAtomWidget wdg(&model);
    auto& ui = wdg.ui;
    auto button = ui->last_picked_element_btn;
    auto button_element_bimap = wdg.m_button_element_bimap;

    std::unordered_set<int> other_button_atomic_numbers;
    for (auto& pair : button_element_bimap) {
        other_button_atomic_numbers.insert(static_cast<int>(pair.right));
    }

    int exp_atomic_number = static_cast<int>(button->getElement());
    for (int atomic_number = 1; atomic_number < 119; ++atomic_number) {
        model.pingValue(ModelKey::ELEMENT, QVariant(atomic_number));
        // The element "last picked" should only correspond to atoms
        // not in the set atom widget panel.
        // No H(1), C(6), N(7), O(8),  F(9), P(15), S(16), Cl(17)
        if (other_button_atomic_numbers.count(atomic_number) == 0) {
            // The element on the "last picked" button should only change if
            // the model's pinged element is not represented by another button
            exp_atomic_number = atomic_number;
        }
        BOOST_TEST(static_cast<int>(button->getElement()) == exp_atomic_number);
        BOOST_TEST(atomic_number_to_name(exp_atomic_number) ==
                   button->toolTip().toStdString());
    }
}

} // namespace sketcher
} // namespace schrodinger
