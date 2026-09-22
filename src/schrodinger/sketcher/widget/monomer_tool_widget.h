#pragma once

#include <memory>
#include <unordered_map>
#include <unordered_set>
#include <string>

#include <boost/bimap.hpp>
#include <boost/signals2/connection.hpp>

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

namespace rdkit_extensions
{
enum class ChainType;
}

namespace sketcher
{

enum class AminoAcidTool;
enum class NucleicAcidTool;
class AminoAcidSymbolPopup;
class CustomNucleotidePopup;
class NucleicAcidSymbolPopup;
class NucleotidePopup;

/**
 * Side bar tool for selecting monomers
 */
class SKETCHER_API MonomerToolWidget : public AbstractDrawToolWidget
{
  public:
    MonomerToolWidget(QWidget* parent = nullptr);
    ~MonomerToolWidget();

    void connectLocalSlots() override;
    void setModel(SketcherModel* model) override;
    void updateCheckedButton() override;
    void updateWidgetsEnabled() override;
    std::unordered_set<QAbstractButton*> getCheckableButtons() override;

  protected:
    /** Rebuild monomer analog popups from the current monomer database. */
    void updateMonomerButtons();

    std::unique_ptr<Ui::MonomerToolWidget> ui;
    boost::bimap<QAbstractButton*, AminoAcidTool> m_button_amino_acid_bimap;
    boost::bimap<QAbstractButton*, NucleicAcidTool> m_button_nucleic_acid_bimap;
    NucleotidePopup* m_rna_popup = nullptr;
    NucleotidePopup* m_dna_popup = nullptr;
    CustomNucleotidePopup* m_custom_nt_popup = nullptr;
    std::unordered_map<QAbstractButton*, AminoAcidSymbolPopup*>
        m_amino_acid_symbol_popups;
    std::unordered_map<QAbstractButton*, NucleicAcidSymbolPopup*>
        m_nucleic_acid_symbol_popups;

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
     * Open the dialog for sketching a custom amino-acid monomer.
     */
    void sketchCustomMonomer();

    /**
     * Respond the the user clicking OK in the custom monomer dialog
     * @param smiles A SMILES string representing the sketched monomer
     * @param monomer_type The monomer type that the user selected in the dialog
     */
    void onCustomMonomerDialogAccepted(
        const std::string& smiles,
        const rdkit_extensions::ChainType monomer_type);

    /**
     * Respond to the user clicking on a specific nucleic acid
     */
    void onNucleicAcidClicked(QAbstractButton* button);

    /**
     * Respond to the user clicking on a specific monomeric connection button
     */
    void onConnectionButtonClicked(int button_id);

  private:
    boost::signals2::scoped_connection m_database_connection;
};

} // namespace sketcher
} // namespace schrodinger
