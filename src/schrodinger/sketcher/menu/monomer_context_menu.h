#pragma once

#include <unordered_set>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/menu/abstract_context_menu.h"

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * Context menu for monomer beads. Currently exposes only Delete; later
 * phases of SKETCH-2648 will add monomer-type-specific actions.
 */
class SKETCHER_API MonomerContextMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    MonomerContextMenu(QWidget* parent = nullptr);

  signals:
    void deleteRequested(const std::unordered_set<const RDKit::Atom*>& atoms);

  private:
    void updateActions() override{};
};

} // namespace sketcher
} // namespace schrodinger
