#pragma once

#include <memory>
#include <unordered_set>

#include <boost/bimap.hpp>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/abstract_draw_tool_widget.h"

class QAbstractButton;
class QWidget;

namespace Ui
{
class MonomerToolWidget;
}

namespace schrodinger
{
namespace sketcher
{

enum class AminoAcidTool;
enum class NucleicAcidTool;
class CustomNucleotidePopup;
class NucleotidePopup;

/**
 * Side bar tool for selecting monomers
 */
class SKETCHER_API MonomerToolWidget : public AbstractDrawToolWidget
{
  public:
    MonomerToolWidget(QWidget* parent = nullptr);
    ~MonomerToolWidget();

    void setModel(SketcherModel* model) override;
    void updateCheckedButton() override;
    std::unordered_set<QAbstractButton*> getCheckableButtons() override;

  protected:
    std::unique_ptr<Ui::MonomerToolWidget> ui;
    boost::bimap<QAbstractButton*, AminoAcidTool> m_button_amino_acid_bimap;
    boost::bimap<QAbstractButton*, NucleicAcidTool> m_button_nucleic_acid_bimap;
    NucleotidePopup* m_rna_popup = nullptr;
    NucleotidePopup* m_dna_popup = nullptr;
    CustomNucleotidePopup* m_custom_nt_popup = nullptr;

    /**
     * Respond to the AMINO or NUCLEIC buttons being clicked, which toggles the
     * tools to the appropriate type of monomer
     */
    void onAminoOrNucleicBtnClicked(QAbstractButton* button);

    /**
     * Respond to the user clicking on a specific amino acid
     */
    void onAminoAcidClicked(QAbstractButton* button);

    /**
     * Respond to the user clicking on a specific nucleic acid
     */
    void onNucleicAcidClicked(QAbstractButton* button);
};

} // namespace sketcher
} // namespace schrodinger
