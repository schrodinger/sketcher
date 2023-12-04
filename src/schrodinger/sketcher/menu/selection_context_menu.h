#pragma once

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace rdkit_extensions
{
enum class Format;
}
namespace sketcher
{

class CutCopyActionManager;
class ModifyAtomsMenu;
class ModifyBondsMenu;
class SketcherModel;
enum class SceneSubset;

class SKETCHER_API SelectionContextMenu : public QMenu
{
    Q_OBJECT
  public:
    SelectionContextMenu(SketcherModel* model, QWidget* parent = nullptr);

    void setConvertToBracketGroupEnabled(bool b);
    void setFlipEnabled(bool b);

    ModifyAtomsMenu* m_modify_atoms_menu = nullptr;
    ModifyBondsMenu* m_modify_bonds_menu = nullptr;

  signals:
    void flipRequested();
    void cutRequested(rdkit_extensions::Format format);
    void copyRequested(rdkit_extensions::Format format, SceneSubset subset);
    void deleteRequested();
    void variableAttachmentBondRequested();
    void bracketSubgroupDialogRequested();
    void invertSelectionRequested();
    void newRGroupRequested();
    void existingRGroupRequested(unsigned int rgroup_number);

  protected:
    void showEvent(QShowEvent* event) override;
    SketcherModel* m_sketcher_model = nullptr;
    CutCopyActionManager* m_cut_copy_actions = nullptr;
    QAction* m_variable_bond_action = nullptr;

  protected slots:
    virtual void updateActionsEnabled();

  private:
    QMenu* createAddToSelectionMenu(SketcherModel* model);
    QMenu* createReplaceSelectionWithMenu(SketcherModel* model);
    QAction* m_bracket_group_action = nullptr;
    QAction* m_flip_action = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
