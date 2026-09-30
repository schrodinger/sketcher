#pragma once

#include <memory>
#include <string>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/modal_dialog.h"

namespace RDKit
{
class ROMol;
}

namespace Ui
{
class CustomMonomerDialog;
}

namespace schrodinger
{

namespace rdkit_extensions
{
enum class ChainType;
}
} // namespace schrodinger

Q_DECLARE_METATYPE(schrodinger::rdkit_extensions::ChainType);

namespace schrodinger
{
namespace sketcher
{

/**
 * Dialog for drawing a custom monomer.
 */
class SKETCHER_API CustomMonomerDialog : public ModalDialog
{
    Q_OBJECT

  public:
    CustomMonomerDialog(const rdkit_extensions::ChainType chain_type,
                        QWidget* parent = nullptr);
    ~CustomMonomerDialog();

    /**
     * Specify the numbered attachment points that must remain in the monomer.
     * If the user removes any of these attachment points and then clicks OK,
     * they will be warned that continuing will remove connections from the
     * monomer.
     */
    void
    setRequiredAttachmentPoints(std::vector<int> required_attachment_points);

    /**
     * Load the specified molecule into the dialog's Sketcher workspace
     */
    void addSMILES(const std::string& smiles);

    /**
     * Overridden QDialog method
     */
    void accept() override;

  signals:
    /**
     * Emitted when the dialog is accepted
     * @param smiles A SMILES string representing the sketched monomer
     * @param monomer_type The monomer type specified when the dialog was opened
     */
    void customMonomerAccepted(const std::string& smiles,
                               const rdkit_extensions::ChainType);

  protected:
    /**
     * Update whether the OK button is enabled. The button is only enabled when
     * the SketcherWidget contains a single molecule that can successfully be
     * converted to a SMILES string.
     */
    void updateOkButton();

    std::unique_ptr<Ui::CustomMonomerDialog> ui;
    rdkit_extensions::ChainType m_chain_type;
    std::vector<int> m_required_attachment_points;
};

} // namespace sketcher
} // namespace schrodinger
