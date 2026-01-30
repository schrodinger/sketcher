#pragma once

#include <unordered_set>

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/menu/abstract_context_menu.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * Context menu for attachment points. Only shows "Delete".
 */
class SKETCHER_API AttachmentPointContextMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    AttachmentPointContextMenu(QWidget* parent = nullptr);

  signals:
    void deleteRequested(const std::unordered_set<const RDKit::Atom*>& atoms,
                         const std::unordered_set<const RDKit::Bond*>& bonds);

  private:
    void updateActions() override{};
};

} // namespace sketcher
} // namespace schrodinger
