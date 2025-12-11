#define BOOST_TEST_MODULE Test_Sketcher

#include "test_abstract_undoable_model.h"

#include <functional>
#include <stdexcept>

#include <QString>
#include <QUndoStack>

#include "../test_common.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

TestUndoableModel::TestUndoableModel(QUndoStack* const undo_stack) :
    AbstractUndoableModel(undo_stack)
{
    connect(this, &TestUndoableModel::inRedo, [this]() { m_in_redo_count++; });
    connect(this, &TestUndoableModel::inUndo, [this]() { m_in_undo_count++; });
    connect(this, &TestUndoableModel::inAdd, [this]() { m_in_add_count++; });
}

void TestUndoableModel::add(int value)
{
    // inAdd won't actually be emitted since AbstractUndoableModel blocks
    // signals outside of redo and undo
    emit inAdd();
    auto redo = [this, value]() {
        m_sum += value;
        emit inRedo();
    };
    auto undo = [this, value]() {
        m_sum -= value;
        emit inUndo();
    };
    QString desc = QString("Add %1").arg(value);
    doCommand(redo, undo, desc);
}

void TestUndoableModel::addTwice(int value)
{
    auto redo = [this, value]() {
        // this will raise an exception since we can't create a command from
        // within a command
        add(value);
        add(value);
    };
    auto undo = [this, value]() {
        add(-value);
        add(-value);
    };
    QString desc = QString("Add %1 twice").arg(value);
    doCommand(redo, undo, desc);
}

void TestUndoableModel::mergeableAdd(int value)
{
    // note that value is passed in as an argument instead of being captured
    auto redo = [this](int value) {
        m_sum += value;
        emit inRedo();
    };
    auto undo = [this](int value) {
        m_sum -= value;
        emit inUndo();
    };
    auto merge = [](int this_value, int other_value) {
        return this_value + other_value;
    };
    QString desc = QString("Mergeable add %1").arg(value);
    int merge_id = 1;
    doMergeableCommand<int>(redo, undo, merge, merge_id, value, desc);
}

BOOST_AUTO_TEST_CASE(test_sumAndSignals)
{
    QUndoStack undo_stack;
    TestUndoableModel model(&undo_stack);
    BOOST_TEST(undo_stack.count() == 0);
    BOOST_TEST(model.m_sum == 0);
    model.add(7);
    BOOST_TEST(undo_stack.count() == 1);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(undo_stack.undoText() == "Add 7");
    BOOST_TEST(model.m_sum == 7);
    BOOST_TEST(model.m_in_add_count == 0);
    BOOST_TEST(model.m_in_redo_count == 1);
    BOOST_TEST(model.m_in_undo_count == 0);
    undo_stack.undo();
    BOOST_TEST(model.m_sum == 0);
    BOOST_TEST(model.m_in_add_count == 0);
    BOOST_TEST(model.m_in_redo_count == 1);
    BOOST_TEST(model.m_in_undo_count == 1);
    undo_stack.redo();
    BOOST_TEST(model.m_sum == 7);
    BOOST_TEST(model.m_in_add_count == 0);
    BOOST_TEST(model.m_in_redo_count == 2);
    BOOST_TEST(model.m_in_undo_count == 1);
    model.add(3);
    BOOST_TEST(undo_stack.count() == 2);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(undo_stack.undoText() == "Add 3");
    BOOST_TEST(model.m_sum == 10);
    BOOST_TEST(model.m_in_add_count == 0);
    BOOST_TEST(model.m_in_redo_count == 3);
    BOOST_TEST(model.m_in_undo_count == 1);
    undo_stack.undo();
    BOOST_TEST(model.m_sum == 7);
    BOOST_TEST(model.m_in_add_count == 0);
    BOOST_TEST(model.m_in_redo_count == 3);
    BOOST_TEST(model.m_in_undo_count == 2);
    undo_stack.undo();
    BOOST_TEST(model.m_sum == 0);
    BOOST_TEST(model.m_in_add_count == 0);
    BOOST_TEST(model.m_in_redo_count == 3);
    BOOST_TEST(model.m_in_undo_count == 3);
    undo_stack.redo();
    BOOST_TEST(model.m_sum == 7);
    BOOST_TEST(model.m_in_add_count == 0);
    BOOST_TEST(model.m_in_redo_count == 4);
    BOOST_TEST(model.m_in_undo_count == 3);
    undo_stack.redo();
    BOOST_TEST(model.m_sum == 10);
    BOOST_TEST(model.m_in_add_count == 0);
    BOOST_TEST(model.m_in_redo_count == 5);
    BOOST_TEST(model.m_in_undo_count == 3);
}

BOOST_AUTO_TEST_CASE(test_createUndoMacro)
{
    QUndoStack undo_stack;
    TestUndoableModel model(&undo_stack);
    BOOST_TEST(undo_stack.count() == 0);
    BOOST_TEST(model.m_sum == 0);
    {
        auto undo_macro = model.createUndoMacro("Macro description");
        model.add(3);
        model.add(5);
        model.add(2);
    }
    BOOST_TEST(undo_stack.count() == 1);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(undo_stack.undoText() == "Macro description");
    BOOST_TEST(model.m_sum == 10);
    undo_stack.undo();
    BOOST_TEST(model.m_sum == 0);
}

/**
 * Test for the move constructor of UndoMacroRAII
 */
BOOST_AUTO_TEST_CASE(test_moveUndoMacro)
{
    QUndoStack undo_stack;
    TestUndoableModel model(&undo_stack);
    BOOST_TEST(undo_stack.count() == 0);
    BOOST_TEST(model.m_sum == 0);
    {
        auto undo_macro = model.createUndoMacro("Macro description");
        auto undo_macro2 = std::move(undo_macro);
        model.add(3);
        model.add(5);
        model.add(2);
    }
    BOOST_TEST(undo_stack.count() == 1);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(undo_stack.undoText() == "Macro description");
    BOOST_TEST(model.m_sum == 10);
    undo_stack.undo();
    BOOST_TEST(model.m_sum == 0);
}

/**
 * Make sure that creating a command from within a command raises an exception
 */
BOOST_AUTO_TEST_CASE(test_exception)
{
#ifdef __APPLE__
    return; // Skip on Mac due to runtime error during exception handling
#endif
    QUndoStack undo_stack;
    TestUndoableModel model(&undo_stack);
    BOOST_CHECK_THROW(model.addTwice(5), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(test_merging)
{
    QUndoStack undo_stack;
    TestUndoableModel model(&undo_stack);
    BOOST_TEST(undo_stack.count() == 0);
    BOOST_TEST(model.m_sum == 0);
    model.mergeableAdd(5);
    BOOST_TEST(undo_stack.count() == 1);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(model.m_sum == 5);
    BOOST_TEST(model.m_in_redo_count == 1);
    BOOST_TEST(model.m_in_undo_count == 0);
    model.mergeableAdd(3);
    BOOST_TEST(undo_stack.count() == 1);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(model.m_sum == 8);
    BOOST_TEST(model.m_in_redo_count == 2);
    BOOST_TEST(model.m_in_undo_count == 0);
    model.mergeableAdd(2);
    BOOST_TEST(undo_stack.count() == 1);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(model.m_sum == 10);
    BOOST_TEST(model.m_in_redo_count == 3);
    BOOST_TEST(model.m_in_undo_count == 0);
    undo_stack.undo();
    BOOST_TEST(model.m_sum == 0);
    BOOST_TEST(model.m_in_redo_count == 3);
    BOOST_TEST(model.m_in_undo_count == 1);
    undo_stack.redo();
    BOOST_TEST(model.m_sum == 10);
    BOOST_TEST(model.m_in_redo_count == 4);
    BOOST_TEST(model.m_in_undo_count == 1);

    model.add(4);
    BOOST_TEST(undo_stack.count() == 2);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(model.m_sum == 14);

    model.mergeableAdd(1);
    BOOST_TEST(undo_stack.count() == 3);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(model.m_sum == 15);
    model.mergeableAdd(7);
    BOOST_TEST(undo_stack.count() == 3);
    BOOST_TEST(undo_stack.canUndo());
    BOOST_TEST(model.m_sum == 22);

    undo_stack.undo();
    BOOST_TEST(model.m_sum == 14);
    undo_stack.undo();
    BOOST_TEST(model.m_sum == 10);
    undo_stack.undo();
    BOOST_TEST(model.m_sum == 0);
    undo_stack.redo();
    BOOST_TEST(model.m_sum == 10);
    undo_stack.redo();
    BOOST_TEST(model.m_sum == 14);
    undo_stack.redo();
    BOOST_TEST(model.m_sum == 22);
}

} // namespace sketcher
} // namespace schrodinger

#include "test_abstract_undoable_model.moc"
