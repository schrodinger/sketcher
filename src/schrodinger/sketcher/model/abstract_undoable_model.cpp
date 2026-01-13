#include "schrodinger/sketcher/model/abstract_undoable_model.h"

#include <stdexcept>

#include <QUndoStack>

#include "schrodinger/sketcher/model/undoable_model_undo_command.h"

namespace schrodinger
{
namespace sketcher
{

AbstractUndoableModel::AbstractUndoableModel(QUndoStack* const undo_stack,
                                             QObject* parent) :
    QObject(parent),
    m_undo_stack(undo_stack)
{
    // prevent signal from being emitted unless we're in a command
    blockSignals(true);
}

void AbstractUndoableModel::setUndoStack(QUndoStack* const undo_stack)
{
    m_undo_stack = undo_stack;
}

void AbstractUndoableModel::doCommand(const std::function<void()> redo,
                                      const std::function<void()> undo,
                                      const QString& description)
{
    throwIfAlreadyAllowingEdits();
    // If error flag was set (e.g., by throwIfAlreadyAllowingEdits), return early
    // to avoid entering Qt code. The error will propagate to the outer doCommand.
    if (m_error_during_command) {
        return;
    }

    // RAII guard to manage command lifetime and ensure proper cleanup
    struct CommandGuard {
        UndoableModelUndoCommand* command;
        bool should_delete = true;
        ~CommandGuard()
        {
            if (should_delete && command) {
                delete command;
            }
        }
    };

    UndoableModelUndoCommand* command =
        new UndoableModelUndoCommand(this, redo, undo, description);
    CommandGuard guard{command};

    m_undo_stack->push(command);

    // Check if error occurred during push (which calls redo/undo callbacks).
    // This allows us to detect errors without throwing through Qt code.
    if (m_error_during_command) {
        QString error = m_error_message;
        m_error_during_command = false;
        m_error_message.clear();

        // Undo the command from the stack. The stack still manages the command's
        // memory, so we shouldn't delete it.
        m_undo_stack->undo();

        // Clear any error flag set during undo
        m_error_during_command = false;
        m_error_message.clear();

        // Stack owns the command now, don't let guard delete it
        guard.should_delete = false;

        throw std::runtime_error(error.toStdString());
    }

    // Success - stack owns the command now, don't let guard delete it
    guard.should_delete = false;
}

void AbstractUndoableModel::throwIfAlreadyAllowingEdits()
{
    if (m_allow_edits) {
        // If we're already allowing edits, then we're either in a command or in
        // the process of creating one.  Creating a command while inside of a
        // command will work correctly during the initial "do" of the command,
        // but will crash when the command is redone.  Set error flags here
        // instead of throwing to avoid exception propagation through Qt code.
        m_error_during_command = true;
        m_error_message = "Cannot create a command while already in edit mode.";
        return;
    }
}

UndoMacroRAII AbstractUndoableModel::createUndoMacro(const QString& description)
{
    return UndoMacroRAII(m_undo_stack, description);
}

void AbstractUndoableModel::beginUndoMacro(const QString& description)
{
    m_undo_stack->beginMacro(description);
}

void AbstractUndoableModel::endUndoMacro()
{
    m_undo_stack->endMacro();
}

UndoMacroRAII::UndoMacroRAII(QUndoStack* undo_stack,
                             const QString& description) :
    m_undo_stack(undo_stack)
{
    m_undo_stack->beginMacro(description);
}

UndoMacroRAII::UndoMacroRAII(UndoMacroRAII&& other)
{
    this->m_undo_stack = other.m_undo_stack;
    other.m_undo_stack = nullptr;
    // The macro has already been begun by other, so we don't want to call
    // beginMacro again
}

UndoMacroRAII::~UndoMacroRAII()
{
    // m_undo_stack will be nullptr if this object has been moved from, in which
    // case the moved-to object is now responsible for ending the macro
    if (m_undo_stack) {
        m_undo_stack->endMacro();
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/model/abstract_undoable_model.moc"