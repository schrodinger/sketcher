#pragma once

#include <unordered_set>
#include <variant>

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
namespace rdkit_extensions
{
enum class Format;
}
namespace sketcher
{

class CutCopyActionManager;
class ExistingRGroupMenu;
class ModifyAtomsMenu;
class ModifyBondsMenu;
class MolModel;
class SketcherModel;
enum class SceneSubset;

class SKETCHER_API SelectionContextMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    SelectionContextMenu(SketcherModel* model, MolModel* mol_model,
                         QWidget* parent = nullptr);

    /**
     * Update the items associated with the context menu and its submenus
     */
    virtual void setContextItems(
        const std::unordered_set<const RDKit::Atom*>& atoms,
        const std::unordered_set<const RDKit::Bond*>& bonds,
        const std::unordered_set<const RDKit::Bond*>& secondary_connections,
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
        const std::unordered_set<const NonMolecularObject*>&
            non_molecular_objects);

    ModifyAtomsMenu* m_modify_atoms_menu = nullptr;
    ModifyBondsMenu* m_modify_bonds_menu = nullptr;

  signals:
    void flipRequested();
    void flipHorizontalRequested();
    void flipVerticalRequested();
    void cutRequested(rdkit_extensions::Format format);
    void copyRequested(rdkit_extensions::Format format, SceneSubset subset);
    void copyAsImageRequested();
    void deleteRequested(
        const std::unordered_set<const RDKit::Atom*>& atoms,
        const std::unordered_set<const RDKit::Bond*>& bonds,
        const std::unordered_set<const RDKit::Bond*>& secondary_connections,
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
        const std::unordered_set<const NonMolecularObject*>&
            non_molecular_objects);
    void variableAttachmentBondRequested(
        const std::unordered_set<const RDKit::Atom*>& atoms);
    void bracketSubgroupDialogRequested(
        const std::unordered_set<const RDKit::Atom*>&);
    void invertSelectionRequested();
    void cleanUpRegionRequested();

  protected:
    SketcherModel* m_sketcher_model = nullptr;
    MolModel* m_mol_model = nullptr;
    CutCopyActionManager* m_cut_copy_actions = nullptr;
    QAction* m_variable_bond_action = nullptr;
    QAction* m_clean_up_region_action = nullptr;

  protected slots:
    virtual void updateActions();

  private:
    QMenu* createAddToSelectionMenu();
    QAction* m_bracket_group_action = nullptr;
    QAction* m_flip_action = nullptr;
    QMenu* m_flip_molecule_menu = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
