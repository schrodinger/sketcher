#pragma once

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"

enum class SceneSubset;

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
    CutCopyActionManager* m_cut_copy_actions = nullptr;
    QAction* m_variable_bond_action = nullptr;

  protected slots:
    void updateActionsEnabled();

  private:
    QMenu* createAddToSelectionMenu(SketcherModel* model);
    QMenu* createReplaceSelectionWithMenu(SketcherModel* model);
    QAction* m_bracket_group_action = nullptr;
    QAction* m_flip_action = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
