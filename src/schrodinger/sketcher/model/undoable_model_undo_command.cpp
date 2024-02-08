#include "schrodinger/sketcher/model/undoable_model_undo_command.h"
#include "schrodinger/sketcher/model/abstract_undoable_model.h"
#include "schrodinger/sketcher/model/tags.h"

#include <QString>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/Geometry/point.h>

#include "schrodinger/sketcher/model/non_molecular_object.h"

namespace schrodinger
{
namespace sketcher
{

template <typename T>
AbstractUndoableModelUndoCommand<T>::AbstractUndoableModelUndoCommand(
    AbstractUndoableModel* const model, const T& redo, const T& undo,
    const QString& description, QUndoCommand* parent) :
    QUndoCommand(description, parent),
    m_model(model),
    m_redo(redo),
    m_undo(undo)
{
}

template <typename T> void AbstractUndoableModelUndoCommand<T>::redo()
{
    do_func(m_redo);
}

template <typename T> void AbstractUndoableModelUndoCommand<T>::undo()
{
    do_func(m_undo);
}

template <typename T>
void AbstractUndoableModelUndoCommand<T>::do_func(const T& func)
{
    bool signals_blocked = m_model->blockSignals(false);
    bool was_allowing_edits = m_model->m_allow_edits;
    m_model->m_allow_edits = true;
    call_func(func);
    m_model->blockSignals(signals_blocked);
    m_model->m_allow_edits = was_allowing_edits;
}

UndoableModelUndoCommand::UndoableModelUndoCommand(
    AbstractUndoableModel* const model, const std::function<void()>& redo,
    const std::function<void()>& undo, const QString& description,
    QUndoCommand* parent) :
    AbstractUndoableModelUndoCommand(model, redo, undo, description, parent)
{
}

void UndoableModelUndoCommand::call_func(const std::function<void()>& func)
{
    func();
}

template <typename T>
UndoableModelMergeableUndoCommand<T>::UndoableModelMergeableUndoCommand(
    AbstractUndoableModel* const model, const std::function<void(T)>& redo,
    const std::function<void(T)>& undo,
    const std::function<T(T, T)>& merge_func, const int merge_id,
    const T init_data, const QString& description, QUndoCommand* parent) :
    AbstractUndoableModelUndoCommand<std::function<void(T)>>(
        model, redo, undo, description, parent),
    m_merge_func(merge_func),
    m_merge_id(merge_id),
    m_data(init_data)
{
}

template <typename T> int UndoableModelMergeableUndoCommand<T>::id() const
{
    return m_merge_id;
}

template <typename T> T UndoableModelMergeableUndoCommand<T>::getData() const
{
    return m_data;
}

template <typename T> bool
UndoableModelMergeableUndoCommand<T>::mergeWith(const QUndoCommand* command)
{
    if (m_merge_id != command->id()) {
        return false;
    }
    auto* mergeable_command =
        dynamic_cast<const UndoableModelMergeableUndoCommand<T>*>(command);
    if (mergeable_command != nullptr) {
        m_data = m_merge_func(m_data, mergeable_command->getData());
    }
    return mergeable_command != nullptr;
}

template <typename T> void UndoableModelMergeableUndoCommand<T>::call_func(
    const std::function<void(T)>& func)
{
    func(m_data);
}

// explicit template instantiations
template class UndoableModelMergeableUndoCommand<int>;
template class UndoableModelMergeableUndoCommand<float>;
template class UndoableModelMergeableUndoCommand<RDGeom::Point3D>;
template class UndoableModelMergeableUndoCommand<
    std::pair<float, RDGeom::Point3D>>;
template class UndoableModelMergeableUndoCommand<
    std::tuple<std::vector<AtomTag>, std::vector<RDGeom::Point3D>,
               std::vector<RDGeom::Point3D>, std::vector<NonMolecularTag>,
               std::vector<RDGeom::Point3D>, std::vector<RDGeom::Point3D>>>;

} // namespace sketcher
} // namespace schrodinger
