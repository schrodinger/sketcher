#pragma once

#include <unordered_set>

#include <QMenu>

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

class ExistingRGroupMenu;
class MolModel;
class ReplaceAtomsWithMenu;
class SetAtomMenuWidget;
class SketcherModel;
enum class AtomQuery;
enum class Element;
enum class ModelKey;

class SKETCHER_API ModifyAtomsMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    ModifyAtomsMenu(SketcherModel* model, MolModel* mol_model,
                    QWidget* parent = nullptr);

  signals:
    /**
     * @param atoms the atoms that should have their charge changed
     * @param increment_by the amount to change the charge by
     */
    void
    adjustChargeRequested(const std::unordered_set<const RDKit::Atom*>& atoms,
                          int increment_by);

    /**
     * @param atoms the atoms that should have their explicit hydrogens either
     * added if there are explicit hydrogens to add, or otherwise removed
     */
    void addRemoveExplicitHydrogensRequested(
        const std::unordered_set<const RDKit::Atom*>& atoms);

    /**
     * @param atoms the atoms that should have their radical electrons changed
     * @param increment_by the amount to change the radical electrons by
     */
    void adjustRadicalElectronsRequested(
        const std::unordered_set<const RDKit::Atom*>& atoms, int increment_by);

    /**
     * @param atoms the atoms that should be changed into rgroups each with the
     * next available rgroup number
     */
    void
    newRGroupRequested(const std::unordered_set<const RDKit::Atom*>& atoms);

    /**
     * @param atoms the atoms that should be changed into rgroups
     * @param rgroup_number the rgroup number to assign when changing the atoms
     */
    void
    existingRGroupRequested(const std::unordered_set<const RDKit::Atom*>& atoms,
                            unsigned int rgroup_number);

    /**
     * @param atoms the atoms that should be changed into a given query type
     * @param query_type the query type to assign to the atoms
     */
    void
    changeTypeRequested(const std::unordered_set<const RDKit::Atom*>& atoms,
                        AtomQuery query_type);

    /**
     * @param atoms the atoms that should have their element changed
     * @param element the element to assign to the atoms
     */
    void
    requestElementChange(const std::unordered_set<const RDKit::Atom*>& atoms,
                         Element element);

    /**
     * Signal to request the scene to show the "Edit Atom Properties" dialog.
     *
     * @param set_to_allowed_list If `true`, show the dialog with the "Allowed
     * List" option set in the "Elements:" combo.
     */
    void showEditAtomPropertiesRequested(const RDKit::Atom* const atom,
                                         const bool set_to_allowed_list);

  protected:
    SketcherModel* m_sketcher_model = nullptr;
    QAction* m_edit_atom_properties_act = nullptr;
    SketcherModel* m_set_atom_model = nullptr;
    ReplaceAtomsWithMenu* m_replace_with_menu = nullptr;
    QAction* m_increase_charge_act = nullptr;
    QAction* m_decrease_charge_act = nullptr;
    QAction* m_add_remove_explicit_h_act = nullptr;
    QAction* m_add_unpaired_electrons_act = nullptr;
    QAction* m_remove_unpaired_electrons_act = nullptr;

  protected:
    /**
     * Disable menu actions that don't make sense given the selected atoms.
     * That is, disable removing unpaired electrons if none are present, etc.
     */
    virtual void updateActions();

    void onSetAtomModelPinged(ModelKey key, QVariant value);

  private:
    QMenu* createElementMenu();
};

class ReplaceAtomsWithMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    ReplaceAtomsWithMenu(MolModel* mol_model, QWidget* parent = nullptr);

    QAction* m_allowed_list_act = nullptr;

  signals:
    void changeTypeRequested(AtomQuery query_type);
    void showEditAtomPropertiesDialogWithAllowedListRequested();
    void newRGroupRequested();
    void existingRGroupRequested(unsigned int rgroup_number);

  protected:
    virtual void updateActions();

  private:
    QMenu* createWildcardMenu();
    ExistingRGroupMenu* m_existing_rgroup_menu = nullptr;
    MolModel* m_mol_model = nullptr;
};

class SKETCHER_API ExistingRGroupMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    ExistingRGroupMenu(MolModel* model, QWidget* parent = nullptr);

  signals:
    void existingRGroupRequested(unsigned int rgroup_number);

  protected:
    virtual void updateActions();

  private:
    MolModel* m_mol_model = nullptr;
};

class AtomContextMenu : public ModifyAtomsMenu
{
    Q_OBJECT
  public:
    AtomContextMenu(SketcherModel* model, MolModel* mol_model,
                    QWidget* parent = nullptr);

  signals:
    void bracketSubgroupDialogRequested(
        const std::unordered_set<const RDKit::Atom*>& atoms);
    void deleteRequested(const std::unordered_set<const RDKit::Atom*>& atoms);

  protected:
    virtual void updateActions();

  private:
    QAction* m_add_brackets_act = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
