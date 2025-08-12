#pragma once

#include "schrodinger/sketcher/model/abstract_undoable_model.h"

class QUndoStack;

namespace schrodinger
{
namespace sketcher
{

class TestUndoableModel : public AbstractUndoableModel
{
    Q_OBJECT
  public:
    TestUndoableModel(QUndoStack* const undo_stack);
    void add(int value);
    void addTwice(int value);
    void mergeableAdd(int value);

    int m_sum = 0;
    int m_in_redo_count = 0;
    int m_in_undo_count = 0;
    int m_in_add_count = 0;
  signals:
    void inRedo();
    void inUndo();
    void inAdd();
};

} // namespace sketcher
} // namespace schrodinger
