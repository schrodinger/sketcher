#pragma once

#include <memory>
#include <unordered_set>

#include "schrodinger/sketcher/dialog/modal_dialog.h"
#include "schrodinger/sketcher/definitions.h"

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
class SKETCHER_API BracketSubgroupDialog : public ModalDialog
{
    Q_OBJECT

  public:
    BracketSubgroupDialog(QWidget* parent = nullptr);
    ~BracketSubgroupDialog();
    void accept() override;
    void setSubgroupType(SubgroupType subgroup_type);
    void setRepeatPattern(RepeatPattern repeat_pattern);
    void setPolymerLabel(const QString& text);
    void setAtoms(std::unordered_set<sketcherAtom*> atoms);

  protected:
    std::unique_ptr<Ui::BracketSubgroupDialog> ui;
    void updateWidgets();
    SubgroupType getSubgroupType() const;
    RepeatPattern getRepeatPattern() const;
    QString getPolymerLabel() const;
    NumericLabelValidator* m_validator = nullptr;

    std::unordered_set<sketcherAtom*> m_atoms;

  signals:
    void bracketSubgroupAccepted(SubgroupType subgroup_type,
                                 RepeatPattern repeat_pattern,
                                 QString polymer_label,
                                 std::unordered_set<sketcherAtom*> atoms);
};

} // namespace sketcher
} // namespace schrodinger
