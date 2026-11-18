#pragma once

#include <memory>
#include <optional>
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
struct MonomerID;
} // namespace rdkit_extensions

namespace sketcher
{

enum class AminoAcidTool;
enum class MonomerToolType;
enum class NucleicAcidTool;
class AminoAcidSymbolPopup;
class CustomNucleotidePopup;
class MonomerDatabasePopup;
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

    /**
     * Activate a database monomer, or use it to mutate the canvas selection.
     * The side-bar page switches to match the monomer: PEPTIDE monomers use
     * the amino acid page, RNA monomers use the nucleic acid page, and CHEM
     * monomers use chem_page. The monomer becomes that page's last-used
     * database monomer.
     *  - Without a canvas selection, the monomer becomes the active draw tool,
     *    and so the page's last monomeric draw tool.
     *  - With a canvas selection, the selected monomers of a matching type are
     *    mutated, and the draw tool is left unchanged.
     * Either way, the page's AMINO_ACID_TOOL or NUCLEIC_ACID_TOOL is left
     * alone, so it still names the page's last standard monomer button.
     *
     * @param monomer The database monomer. Nothing happens if its symbol is
     * empty.
     * @param chem_page The page to use for a CHEM monomer, which belongs on
     * either page. Defaults to the current page (MONOMER_TOOL_TYPE). Pass it
     * when restoring a page's monomer before the model has switched to that
     * page.
     */
    void activateMonomerDbMonomer(
        const rdkit_extensions::MonomerID& monomer,
        std::optional<MonomerToolType> chem_page = std::nullopt);

    /**
     * Set whether the Monomer Database popup's Load Database... (or Update
     * Database...) link is visible. This applies even if the popup hasn't
     * been created yet.
     */
    void setLoadMonomerDatabaseVisible(bool visible);

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
     * Respond to the user clicking the Monomer Database button by showing the
     * Monomer Database popup, or by closing it if it's already open
     */
    void onMonomerDbButtonClicked();

    /**
     * Respond to the user clicking the last-used database monomer button by
     * activating the current page's last-used database monomer, if there is
     * one
     */
    void onMonomerDbLastUsedButtonClicked();

    /**
     * @return The Monomer Database popup, which is created the first time
     * this is called
     */
    MonomerDatabasePopup* getMonomerDatabasePopup();

    /**
     * Respond to the user clicking on a specific monomeric connection button
     */
    void onConnectionButtonClicked(int button_id);

  private:
    boost::signals2::scoped_connection m_database_connection;
    MonomerDatabasePopup* m_monomer_db_popup = nullptr;
    bool m_load_monomer_database_visible = true;
};

} // namespace sketcher
} // namespace schrodinger
