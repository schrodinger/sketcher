#include "schrodinger/sketcher/menu/abstract_context_menu.h"

#include "schrodinger/sketcher/model/non_molecular_object.h"
namespace schrodinger
{
namespace sketcher
{

AbstractContextMenu::AbstractContextMenu(QWidget* parent) : QMenu(parent)
{
}

void AbstractContextMenu::setContextItems(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const std::unordered_set<const RDKit::Bond*>& secondary_connections,
    const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{
    m_atoms = atoms;
    m_bonds = bonds;
    m_secondary_connections = secondary_connections;
    m_sgroups = sgroups;
    m_non_molecular_objects = non_molecular_objects;
}

void AbstractContextMenu::showEvent(QShowEvent* event)
{
    updateActions();
    QMenu::showEvent(event);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/abstract_context_menu.moc"