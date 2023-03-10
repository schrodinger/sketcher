#include "schrodinger/sketcher/molviewer/undoable_model_undo_command.h"

#include <QString>

namespace schrodinger
{
namespace sketcher
{

UndoableModelUndoCommand::UndoableModelUndoCommand(
    AbstractUndoableModel* const model, const std::function<void()>& redo,
    const std::function<void()>& undo, const QString& description,
    QUndoCommand* parent) :
    QUndoCommand(description, parent),
    m_model(model),
    m_redo(redo),
    m_undo(undo)
{
}

void UndoableModelUndoCommand::redo()
{
    do_func(m_redo);
}

void UndoableModelUndoCommand::undo()
{
    do_func(m_undo);
}
void UndoableModelUndoCommand::do_func(const std::function<void()>& func)
{
    bool signals_blocked = m_model->blockSignals(false);
    bool was_in_command = m_model->m_in_command;
    m_model->m_in_command = true;
    func();
    m_model->blockSignals(signals_blocked);
    m_model->m_in_command = was_in_command;
}

} // namespace sketcher
} // namespace schrodinger
