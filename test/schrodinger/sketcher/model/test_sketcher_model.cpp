#define BOOST_TEST_MODULE Test_Sketcher
#include <memory>

#include <QGraphicsTextItem>
#include <QSignalSpy>
#include <boost/test/unit_test.hpp>

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/reaction_arrow.h"
#include "test_common.h"
#include "test_sketcherScene.h"

Q_DECLARE_METATYPE(schrodinger::sketcher::ModelKey);
Q_DECLARE_METATYPE(std::unordered_set<schrodinger::sketcher::ModelKey>);
Q_DECLARE_METATYPE(std::unordered_set<QGraphicsItem*>);

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

BOOST_TEST_DONT_PRINT_LOG_VALUE(QVariant)

namespace schrodinger
{
namespace sketcher
{

/**
 * Subclass that exposes protected data for testing purposes.
 */
class TestSketcherModel : public SketcherModel
{
  public:
    std::unordered_map<ModelKey, QVariant>& getModelMap()
    {
        return m_model_map;
    }
};

/**
 * Verify that the return values of `get_model_keys()` match the keys in the
 * model dictionary.
 */
BOOST_AUTO_TEST_CASE(model_keys)
{

    TestSketcherModel model;
    auto model_map = model.getModelMap();
    auto all_keys = get_model_keys();
    BOOST_TEST(model_map.size() == all_keys.size());
    for (auto& key : all_keys) {
        bool element_found = true;
        try {
            model_map.at(key);
        } catch (const std::out_of_range&) {
            element_found = false;
        }
        BOOST_TEST(element_found);
    }
}

/**
 * Verify correct behavior of getter method, setter method, and signals.
 */
BOOST_AUTO_TEST_CASE(get_set_signal)
{

    qRegisterMetaType<ModelKey>("ModelKey");
    qRegisterMetaType<std::unordered_set<ModelKey>>(
        "std::unordered_set<ModelKey>");
    auto var = QVariant::fromValue(ModelKey::ATOM_TOOL);
    BOOST_TEST(var.typeId() >= QMetaType::User);
    BOOST_TEST(static_cast<int>(var.value<ModelKey>()) ==
               static_cast<int>(ModelKey::ATOM_TOOL));

    std::unordered_map<ModelKey, QVariant> new_values = {
        {ModelKey::NEW_STRUCTURES_REPLACE_CONTENT, false},
        {ModelKey::SHOW_LID_LEGEND, true},
        {ModelKey::ALLOW_MULTIPLE_RXNS, true},
        {ModelKey::SHOW_VALENCE_ERRORS, false},
        {ModelKey::COLOR_HETEROATOMS, false},
        {ModelKey::SHOW_STEREOCENTER_LABELS, false},
        {ModelKey::USE_IMPLICIT_HYDROGENS, true},
        {ModelKey::SELECTION_TOOL, QVariant::fromValue(SelectionTool::LASSO)},
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MOVE_ROTATE)},
        {ModelKey::ATOM_TOOL, QVariant::fromValue(AtomTool::QUERY)},
        {ModelKey::BOND_TOOL, QVariant::fromValue(BondTool::DOUBLE)},
        {ModelKey::CHARGE_TOOL, QVariant::fromValue(ChargeTool::INCREASE)},
        {ModelKey::RING_TOOL, QVariant::fromValue(RingTool::CYCLOBUTANE)},
        {ModelKey::ENUMERATION_TOOL,
         QVariant::fromValue(EnumerationTool::EXISTING_RGROUP)},
        {ModelKey::ELEMENT, QVariant::fromValue(Element::N)},
        {ModelKey::ATOM_QUERY, QVariant::fromValue(AtomQuery::Q)},
        {ModelKey::RGROUP_NUMBER, 2u},
        {ModelKey::RESIDUE_TYPE, QString("foo")},
    };

    TestSketcherModel model;
    QSignalSpy pinged_spy(&model, &SketcherModel::valuePinged);
    QSignalSpy changed_spy(&model, &SketcherModel::valuesChanged);
    std::vector<QSignalSpy*> spies = {&pinged_spy, &changed_spy};
    for (auto& key : get_model_keys()) {
        if (key == ModelKey::RESIDUE_TYPE) {
            // This value is not stored as an int-like object
            continue;
        }
        // Access the current value
        const auto& value = model.getValue(key);
        for (auto spy : spies) {
            BOOST_TEST(spy->empty());
        }
        BOOST_TEST(value.isValid());

        // Assign a new value
        auto new_value = new_values[key];
        model.setValue(key, new_value);

        // Verify that a signal was emitted
        for (auto spy : spies) {
            BOOST_TEST(spy->count() == 1);
        }

        {
            // valuePinged() is emitted with a key and corresponding value
            auto args = pinged_spy.takeLast();
            BOOST_TEST(args.at(0).typeId() >= QMetaType::User);
            auto signal_key = args.at(0).value<ModelKey>();
            BOOST_TEST(static_cast<int>(signal_key) == static_cast<int>(key));
            BOOST_TEST(args.at(1).toInt() == new_value.toInt());
        }

        {
            // valuesChanged() is emitted with a set of changed keys
            auto args = changed_spy.takeLast();
            BOOST_TEST(args.at(0).typeId() >= QMetaType::User);
            auto keys = args.at(0).value<std::unordered_set<ModelKey>>();
            BOOST_TEST(keys.size() == 1);
        }

        // Assign the same value again, and verify that the "pinged" signal is
        // emitted but the "changed" signal is not
        model.setValue(key, new_value);
        BOOST_TEST(pinged_spy.count() == 1);
        auto args = pinged_spy.takeLast();
        BOOST_TEST(args.at(0).typeId() >= QMetaType::User);
        auto signal_key = args.at(0).value<ModelKey>();
        BOOST_TEST(static_cast<int>(signal_key) == static_cast<int>(key));
        BOOST_TEST(args.at(1).toInt() == new_value.toInt());
        BOOST_TEST(changed_spy.empty());
    }
}

/**
 * Verify correct behavior of multiple setter method `setValues()`.
 */
BOOST_AUTO_TEST_CASE(setValues)
{

    qRegisterMetaType<std::unordered_set<ModelKey>>(
        "std::unordered_set<ModelKey>");
    TestSketcherModel model;
    QSignalSpy changed_spy(&model, &SketcherModel::valuesChanged);

    std::unordered_set<ModelKey> keys_to_change = {ModelKey::SELECTION_TOOL,
                                                   ModelKey::DRAW_TOOL};
    std::unordered_set<ModelKey> keys_stay_the_same = {
        ModelKey::ENUMERATION_TOOL, ModelKey::ELEMENT};
    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::SELECTION_TOOL, QVariant::fromValue(SelectionTool::LASSO)},
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::BOND)},
    };
    for (auto key : keys_to_change) {
        BOOST_TEST(model.getValue(key) != kv_pairs[key]);
    }
    for (auto key : keys_stay_the_same) {
        kv_pairs.emplace(key, model.getValue(key));
    }

    BOOST_TEST(kv_pairs.size() ==
               keys_to_change.size() + keys_stay_the_same.size());

    model.setValues(kv_pairs);
    BOOST_TEST(changed_spy.count() == 1);
    auto args = changed_spy.takeLast();
    BOOST_TEST(args.at(0).typeId() >= QMetaType::User);
    auto signal_keys = args.at(0).value<std::unordered_set<ModelKey>>();

    std::function<bool(ModelKey)> is_in_signal_keys =
        [signal_keys](ModelKey key) { return signal_keys.count(key) == 1; };
    BOOST_TEST(std::all_of(keys_to_change.begin(), keys_to_change.end(),
                           is_in_signal_keys));
    BOOST_TEST(std::none_of(keys_stay_the_same.begin(),
                            keys_stay_the_same.end(), is_in_signal_keys));

    for (auto& [key, value] : kv_pairs) {
        BOOST_TEST(model.getValue(key) == value);
    }
}

/**
 * Verify correct behavior when interacting with reaction arrow coordinates.
 */
BOOST_AUTO_TEST_CASE(reactions)
{
    testSketcherScene scene;
    auto model = scene.getModel();

    BOOST_TEST(!model->hasReaction());

    QPointF pos(0.0, 0.0);
    scene.addReactionArrowAt(pos);
    BOOST_TEST(model->hasReaction());

    scene.undo();
    BOOST_TEST(!model->hasReaction());
}

/**
 * Verify that the model can retrieve information from the scene undo stack.
 */
BOOST_AUTO_TEST_CASE(getUndoStackData)
{
    testSketcherScene scene;
    auto model = scene.getModel();

    auto [can_undo, can_redo] = model->getUndoStackData();
    BOOST_TEST(!can_undo);
    BOOST_TEST(!can_redo);

    // Do something in the scene so that it gets put onto the undo stack
    scene.importText("CCC");
    std::tie(can_undo, can_redo) = model->getUndoStackData();
    BOOST_TEST(can_undo);
    BOOST_TEST(!can_redo);

    // Undo will move an item from the undo stack to the redo stack
    scene.undo();
    std::tie(can_undo, can_redo) = model->getUndoStackData();
    BOOST_TEST(!can_undo);
    BOOST_TEST(can_redo);

    // Redo will do the opposite
    scene.redo();
    std::tie(can_undo, can_redo) = model->getUndoStackData();
    BOOST_TEST(can_undo);
    BOOST_TEST(!can_redo);
}

} // namespace sketcher
} // namespace schrodinger
