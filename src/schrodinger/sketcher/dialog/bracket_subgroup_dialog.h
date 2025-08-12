#pragma once

#include <memory>
#include <unordered_set>
#include <variant>

#include <rdkit/GraphMol/SubstanceGroup.h>
#include <rdkit/GraphMol/Atom.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/modal_dialog.h"
#include "schrodinger/sketcher/model/mol_model.h"

class sketcherAtom;
namespace Ui
{
class BracketSubgroupDialog;
}

namespace schrodinger
{
namespace sketcher
{

class NumericLabelValidator;
enum class SubgroupType;
enum class RepeatPattern;

/**
 * Dialog for specifying bracket subgroups.
 */
class SKETCHER_API AbstractBracketSubgroupDialog : public ModalDialog
{
    Q_OBJECT

  public:
    AbstractBracketSubgroupDialog(QWidget* parent = nullptr);
    ~AbstractBracketSubgroupDialog();
    void setSubgroupType(SubgroupType subgroup_type);
    void setRepeatPattern(RepeatPattern repeat_pattern);
    void setPolymerLabel(const QString& text);

  protected:
    std::unique_ptr<Ui::BracketSubgroupDialog> ui;
    void updateWidgets();
    SubgroupType getSubgroupType() const;
    RepeatPattern getRepeatPattern() const;
    QString getPolymerLabel() const;
    NumericLabelValidator* m_validator = nullptr;
};

/**
 * Dialog for specifying bracket subgroups when using sketcherScene
 */
class SKETCHER_API BracketSubgroupDialogDeprecated
    : public AbstractBracketSubgroupDialog
{
    Q_OBJECT

  public:
    BracketSubgroupDialogDeprecated(QWidget* parent = nullptr);
    void accept() override;
    void setAtoms(const std::unordered_set<sketcherAtom*>& atoms);

  protected:
    std::unordered_set<sketcherAtom*> m_atoms;

  signals:
    void bracketSubgroupAccepted(SubgroupType subgroup_type,
                                 RepeatPattern repeat_pattern,
                                 QString polymer_label,
                                 std::unordered_set<sketcherAtom*> atoms);
};

/**
 * Dialog for specifying bracket subgroups when using the molviewer Scene
 */
class SKETCHER_API BracketSubgroupDialog : public AbstractBracketSubgroupDialog
{
    Q_OBJECT

  public:
    BracketSubgroupDialog(MolModel* const mol_model, QWidget* parent = nullptr);
    void accept() override;

    /**
     * Specify the atoms that this dialog should affect.  If the dialog is
     * accepted, a new S-group will be created for these atoms.
     */
    void setAtoms(const std::unordered_set<const RDKit::Atom*>& atoms);

    /**
     * Specify the S-group that this dialog should affect.  If the dialog is
     * accepted, this S-group will be modified.
     */
    void setSubgroup(const RDKit::SubstanceGroup* const s_group);

  protected:
    MolModel* m_mol_model;
    std::variant<std::unordered_set<const RDKit::Atom*>,
                 const RDKit::SubstanceGroup*>
        m_atoms_or_s_group;
};

} // namespace sketcher
} // namespace schrodinger
