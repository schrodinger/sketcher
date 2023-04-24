/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <functional>

#include <QObject>
#include <QUndoStack>

#include "schrodinger/sketcher/definitions.h"

class QString;
class QUndoStack;

namespace schrodinger
{
namespace sketcher
{

class UndoableModelUndoCommand;

/**
 * An RAII class for creating undo macros in a QUndoStack
 */
class SKETCHER_API UndoMacroRAII
{
  public:
    /**
     * @param undo_stack The undo stack to create a macro in
     * @param description A description of the macro
     */
    UndoMacroRAII(QUndoStack* undo_stack, const QString& description);

    // move constructor, which will prevent the moved-from object from ending
    // the macro when it's destroyed
    UndoMacroRAII(UndoMacroRAII&& other);

    // make sure we don't inadvertently create a copy of this object
    UndoMacroRAII(const UndoMacroRAII&) = delete;
    UndoMacroRAII& operator=(const UndoMacroRAII&) = delete;

    // don't allow move assignment since we don't have a default constructor, so
    // we can't create an empty object to move into
    UndoMacroRAII& operator=(UndoMacroRAII&& rhs) = delete;

    // the destructor will end the macro, unless this object has been moved from
    ~UndoMacroRAII();

  protected:
    QUndoStack* m_undo_stack;
};

/**
 * A base class for models where all changes are undoable using a QUndoStack.
 * Undoable methods can be defined by calling doCommand with redo and undo
 * lambdas.  For example:
 *
 * void MyUndoableModel::add(int value)
 * {
 *     auto redo = [this, value]()
 *     {
 *         m_sum += value;
 *         emit valueChanged(m_sum);
 *     };
 *     auto undo = [this, value]()
 *     {
 *         m_sum -= value;
 *         emit valueChanged(m_sum);
 *     };
 *     QString desc = QString("Add %1").arg(value);
 *     doCommand(redo, undo, desc);
 * }
 *
 * With the above definitions, any calls to `add` will immediately add `value`
 * to `m_sum` and will be undoable.
 *
 * Note that signals may only be emitted from inside a redo or undo.  Signals
 * emitted at other times will be blocked.  This is to make it easier to catch
 * bugs caused by accidentally emitting signals when a command is created
 * instead of when a command is executed (i.e. accidentally putting the emit
 * outside of the lambda instead of inside).
 *
 * Also note that doCommand must not be called from inside a redo or undo, as
 * this will raise an exception.  Without the exception, this would instead
 * result in a downstream crash.
 *
 */
class SKETCHER_API AbstractUndoableModel : public QObject
{
    Q_OBJECT
  public:
    /**
     * @param undo_stack The undo stack to use.  If not given here, setUndoStack
     * must be called before this instance creates any undo commands.
     */
    AbstractUndoableModel(QUndoStack* undo_stack = nullptr,
                          QObject* parent = nullptr);

    /**
     * Specify the undo stack to use
     */
    void setUndoStack(QUndoStack* const undo_stack);

    /**
     * Return an RAII object that will manage the lifetime of an undo macro
     * @param description A description of the macro
     */
    UndoMacroRAII createUndoMacro(const QString& description);

  protected:
    /**
     * Create an undo command and add it to the undo stack, which automatically
     * executes the command
     *
     * @param redo The function to call on redo (or for the initial do)
     * @param undo The function to call on undo
     * @param description A description of the command
     */
    void doCommand(const std::function<void()> redo,
                   const std::function<void()> undo,
                   const QString& description);

    QUndoStack* m_undo_stack;

    // m_in_command is updated by UndoableModelUndoCommand immediately before
    // and after running a command
    bool m_in_command = false;
    friend class UndoableModelUndoCommand;
};

} // namespace sketcher
} // namespace schrodinger
