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
    if (m_error_during_command) {
        return;  // Nested call detected, propagate error to outer doCommand
    }

    // RAII guard ensures command cleanup on exceptions
    struct CommandGuard {
        UndoableModelUndoCommand* cmd;
        bool should_delete = true;
        ~CommandGuard() { if (should_delete && cmd) delete cmd; }
    };

    UndoableModelUndoCommand* command =
        new UndoableModelUndoCommand(this, redo, undo, description);
    CommandGuard guard{command};

    m_undo_stack->push(command);

    if (m_error_during_command) {
        QString error = m_error_message;
        m_error_during_command = false;
        m_error_message.clear();

        m_undo_stack->undo();  // Remove failed command from stack
        guard.should_delete = false;  // Stack manages command memory

        throw std::runtime_error(error.toStdString());
    }

    guard.should_delete = false;  // Success - stack owns command
}

void AbstractUndoableModel::throwIfAlreadyAllowingEdits()
{
    if (m_allow_edits) {
        // Nested doCommand() detected. Set error flags instead of throwing to
        // avoid exception propagation through Qt's static code on macOS.
        m_error_during_command = true;
        m_error_message = "Cannot create a command while already in edit mode.";
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