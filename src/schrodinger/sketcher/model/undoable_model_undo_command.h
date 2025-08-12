#pragma once

#include <functional>

#include <QUndoCommand>

#include "schrodinger/sketcher/definitions.h"

class QObject;

namespace schrodinger
{
namespace sketcher
{

class AbstractUndoableModel;

/**
 * A base class for code shared by UndoableModelUndoCommand and
 * UndoableModelMergeableUndoCommand
 */
template <typename T> class SKETCHER_API AbstractUndoableModelUndoCommand
    : public QUndoCommand
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
    AbstractUndoableModelUndoCommand(AbstractUndoableModel* const model,
                                     const T& redo, const T& undo,
                                     const QString& description,
                                     QUndoCommand* parent = nullptr);

    // overridden QUndoCommand methods
    void redo() override;
    void undo() override;

  protected:
    /**
     * Run the specified function.  While the function is being run, signals
     * will be unblocked on the associated AbstractUndoableModel instance.
     * AbstractUndoableModel::m_allow_edits will also be updated so that
     * AbstractUndoableModel can catch commands that create commands (which
     * isn't allowed and would lead to a crash).
     *
     * @param func The function to run.  Should be either m_redo or m_undo.
     */
    void do_func(const T& func);

    /**
     * Run the specified function, passing in any necessary arguments
     *
     * @param func The function to run.  Should be either m_redo or m_undo.
     */
    virtual void call_func(const T& func) = 0;

    AbstractUndoableModel* const m_model;
    T m_redo;
    T m_undo;
};

/**
 * An undo command for use with an AbstractUndoableModel.  This class accepts
 * lambdas for redo and undo methods so that we don't need to create a separate
 * undo command class for each undoable method.
 */
class SKETCHER_API UndoableModelUndoCommand
    : public AbstractUndoableModelUndoCommand<std::function<void()>>
{
  public:
    /**
     * See AbstractUndoableModelUndoCommand for parameter documentation
     */
    UndoableModelUndoCommand(AbstractUndoableModel* const model,
                             const std::function<void()>& redo,
                             const std::function<void()>& undo,
                             const QString& description,
                             QUndoCommand* parent = nullptr);

  protected:
    // overridden AbstractUndoableModelUndoCommand method
    void call_func(const std::function<void()>& func) override;
};

/**
 * A mergeable undo command for use with an AbstractUndoableModel.  In addition
 * to the redo and undo lambdas, this class also accepts a lambda for merging
 * commands with the same merge_id.
 */
template <typename T> class SKETCHER_API UndoableModelMergeableUndoCommand
    : public AbstractUndoableModelUndoCommand<std::function<void(T)>>
{
  public:
    /**
     * @param model The undoable model instance that this command will be part
     * of
     * @param redo The redo function that accepts a data argument of type T
     * @param undo The undo function that accepts a data argument of type T
     * @param merge_func A function that accepts two data arguments of type T
     * and returns the merged value
     * @param merge_id The merge ID for the command.  Two commands will only be
     * merged if they have the same merge ID.
     * @param init_data The initial value for the data argument when calling
     * redo and undo.  This value may later be modified as a result of merging
     * commands.
     * @param description A description of this command
     * @param parent The parent command, if any
     */
    UndoableModelMergeableUndoCommand(AbstractUndoableModel* const model,
                                      const std::function<void(T)>& redo,
                                      const std::function<void(T)>& undo,
                                      const std::function<T(T, T)>& merge_func,
                                      const int merge_id, const T init_data,
                                      const QString& description,
                                      QUndoCommand* parent = nullptr);

    /**
     * @return the data that should be passed into the redo, undo, and merge
     * functions
     */
    T getData() const;

    // overridden QUndoCommand methods
    int id() const override;
    bool mergeWith(const QUndoCommand* command) override;

  protected:
    std::function<T(T, T)> m_merge_func;
    int m_merge_id;
    T m_data;

    // overridden AbstractUndoableModelUndoCommand method
    void call_func(const std::function<void(T)>& func) override;
};

} // namespace sketcher
} // namespace schrodinger
