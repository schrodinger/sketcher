#pragma once

#include <unordered_set>

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"

namespace RDKit
{
class Atom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

class ExistingRGroupMenu;
class ReplaceAtomsWithMenu;
class SetAtomMenuWidget;
class SketcherModel;
enum class AtomQuery;
enum class Element;
enum class ModelKey;

class SKETCHER_API ModifyAtomsMenu : public QMenu
{
    Q_OBJECT
  public:
    ModifyAtomsMenu(SketcherModel* model, QWidget* parent = nullptr);

    /**
     * Update the atoms associated with the context menu actions, updating
     * enabled actions as appropriate.
     */
    void setContextAtoms(const std::unordered_set<const RDKit::Atom*>& atoms);

    /**
     * Overrides the show event to also call updateActionsEnabled
     */
    void showEvent(QShowEvent* event);

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
    void showEditAtomPropertiesRequested(bool set_to_allowed_list);

  protected:
    SketcherModel* m_sketcher_model = nullptr;
    QAction* m_edit_atom_properties_act = nullptr;
    SketcherModel* m_set_atom_model = nullptr;
    SetAtomMenuWidget* m_element_widget = nullptr;
    ReplaceAtomsWithMenu* m_replace_with_menu = nullptr;
    QAction* m_increase_charge_act = nullptr;
    QAction* m_decrease_charge_act = nullptr;
    QAction* m_add_remove_explicit_h_act = nullptr;
    QAction* m_add_unpaired_electrons_act = nullptr;
    QAction* m_remove_unpaired_electrons_act = nullptr;
    std::unordered_set<const RDKit::Atom*> m_atoms;

  protected slots:
    /**
     * Disable menu actions that don't make sense given the selected atoms.
     * That is, disable removing unpaired electrons if none are present, etc.
     */
    virtual void updateActionsEnabled();

    void onSetAtomModelPinged(ModelKey key, QVariant value);

  private:
    QMenu* createElementMenu();
};

class ReplaceAtomsWithMenu : public QMenu
{
    Q_OBJECT
  public:
    ReplaceAtomsWithMenu(SketcherModel* model, QWidget* parent = nullptr);
    void showEvent(QShowEvent* event);
    QAction* m_allowed_list_act = nullptr;

  signals:
    void changeTypeRequested(AtomQuery query_type);
    void showEditAtomPropertiesDialogWithAllowedListRequested();
    void newRGroupRequested();
    void existingRGroupRequested(unsigned int rgroup_number);

  private:
    QMenu* createWildcardMenu();
    ExistingRGroupMenu* m_existing_rgroup_menu = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
};

class SKETCHER_API ExistingRGroupMenu : public QMenu
{
    Q_OBJECT
  public:
    ExistingRGroupMenu(SketcherModel* model, QWidget* parent = nullptr);
    void showEvent(QShowEvent* event);

  signals:
    void existingRGroupRequested(unsigned int rgroup_number);

  private:
    SketcherModel* m_sketcher_model = nullptr;
    void updateItems();
};

class AtomContextMenu : public ModifyAtomsMenu
{
    Q_OBJECT
  public:
    AtomContextMenu(SketcherModel* model, QWidget* parent = nullptr);

  signals:
    void bracketSubgroupDialogRequested();
    void deleteRequested(const std::unordered_set<const RDKit::Atom*>& atoms);

  protected slots:
    void updateActionsEnabled();

  private:
    QAction* m_add_brackets_act = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
