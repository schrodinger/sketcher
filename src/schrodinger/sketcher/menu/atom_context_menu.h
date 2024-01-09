#pragma once

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"

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

    ReplaceAtomsWithMenu* m_replace_with_menu = nullptr;
    void showEvent(QShowEvent* event);

    /**
     * Disable menu actions that don't make sense given the selected atoms.
     * That is, disable removing unpaired electrons if none are present, etc.
     */
    virtual void updateActionsEnabled();

  signals:
    void increaseChargeRequested();
    void decreaseChargeRequested();
    void addRemoveExplicitHydrogensRequested();
    void addUnpairedElectronRequested();
    void removeUnpairedElectronRequested();
    void newRGroupRequested();
    void existingRGroupRequested(unsigned int rgroup_number);

    /**
     * Signal to request the scene to change the element for all atoms
     * associated with the context menu.
     *
     * @param element The element to assign to the atoms
     */
    void requestContextMenuElementChange(Element element);

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
    QAction* m_increase_charge_act = nullptr;
    QAction* m_decrease_charge_act = nullptr;
    QAction* m_add_remove_explicit_h_act = nullptr;
    QAction* m_add_unpaired_electrons_act = nullptr;
    QAction* m_remove_unpaired_electrons_act = nullptr;

  protected slots:
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
    void setAddToBracketGroupEnabled(bool enabled);

  signals:
    void deleteRequested();
    void existingRGroupRequested(unsigned int rgroup_number);
    void bracketSubgroupDialogRequested();

  private:
    QAction* m_add_brackets_act = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
