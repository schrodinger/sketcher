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
    throwIfInCommand();
    UndoableModelUndoCommand* command =
        new UndoableModelUndoCommand(this, redo, undo, description);
    m_undo_stack->push(command);
}

void AbstractUndoableModel::throwIfInCommand()
{
    if (m_in_command) {
        // Creating a command while inside of a command will work correctly
        // during the initial "do" of the command, but will crash when the
        // command is redone.  Throw an exception here so the problem is
        // immediately obvious.
        throw std::runtime_error(
            "Cannot create a command from within a command");
    }
}

UndoMacroRAII AbstractUndoableModel::createUndoMacro(const QString& description)
{
    return UndoMacroRAII(m_undo_stack, description);
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