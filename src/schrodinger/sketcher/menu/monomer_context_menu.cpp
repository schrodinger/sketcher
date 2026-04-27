#include "schrodinger/sketcher/menu/monomer_context_menu.h"

namespace schrodinger
{
namespace sketcher
{

MonomerContextMenu::MonomerContextMenu(QWidget* parent) :
    AbstractContextMenu(parent)
{
    setTitle("Monomer");
    addAction("Delete", this, [this]() { emit deleteRequested(m_atoms); });
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/monomer_context_menu.moc"
