#include "schrodinger/sketcher/menu/attachment_point_context_menu.h"

namespace schrodinger
{
namespace sketcher
{

AttachmentPointContextMenu::AttachmentPointContextMenu(QWidget* parent) :
    AbstractContextMenu(parent)
{
    setTitle("Attachment Point");
    addAction("Delete", this,
              [this]() { emit deleteRequested(m_atoms, m_bonds); });
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/attachment_point_context_menu.moc"
