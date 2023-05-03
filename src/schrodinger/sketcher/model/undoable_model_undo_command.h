#pragma once

#include <functional>

#include <QUndoCommand>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/abstract_undoable_model.h"

class QObject;

namespace schrodinger
{
namespace sketcher
{

/**
 * An undo command for use with an AbstractUndoableModel.  This class accepts
 * lambdas for redo and undo methods so that we don't need to create a separate
 * undo command class for each undoable method.
 */
class SKETCHER_API UndoableModelUndoCommand : public QUndoCommand
{
  public:
    /**
     * @param model The undoable model instance that this command will be part
     * of
     * @param redo The redo function
     * @param undo The undo function
     * @param description A description of this command
     * @param parent The parent command, if any
     */
    UndoableModelUndoCommand(AbstractUndoableModel* const model,
                             const std::function<void()>& redo,
                             const std::function<void()>& undo,
                             const QString& description,
                             QUndoCommand* parent = nullptr);

    void redo() override;
    void undo() override;

  protected:
    /**
     * Run the specified function.  While the function is being run, signals
     * will be unblocked on the associated AbstractUndoableModel instance.
     * AbstractUndoableModel::m_in_command will also be updated so that
     * AbstractUndoableModel can catch commands that create commands (which
     * isn't allowed and would lead to a crash).
     *
     * @param func The function to run.  Should be either m_redo or m_undo.
     */
    void do_func(const std::function<void()>& func);

    AbstractUndoableModel* const m_model;
    std::function<void()> m_redo;
    std::function<void()> m_undo;
};

} // namespace sketcher
} // namespace schrodinger
