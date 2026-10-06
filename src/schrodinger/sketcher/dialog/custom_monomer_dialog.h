#pragma once

#include <memory>
#include <string>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/resizable_model_dialog.h"

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
 *
 * To ensure that this dialog can fit into Live Design's Sketcher frame, its
 * minimum height must be no larger than a standard SketcherWidget's minimum
 * height.  To accomplish this, we hide the atomistic/monomeric toggle and
 * shrink the typical ModalDialog border. Additionally, if the dialog is small
 * enough, we allow the title bar to be shrunk and allow the bottom bar to be
 * moved into the embedded SketcherWidget underneath the View (i.e. next to the
 * bottom of the side bar instead of below the side bar).
 */
class SKETCHER_API CustomMonomerDialog : public ResizableModelDialog
{
    Q_OBJECT

  public:
    CustomMonomerDialog(const rdkit_extensions::ChainType chain_type,
                        QWidget* parent = nullptr);
    ~CustomMonomerDialog();

    /**
     * Overridden Qt methods to let Qt know that the dialog can be made smaller
     * by moving the bottom bar into the embedded SketcherWidget, but that the
     * dialog's preferred size leaves enough room to keep the bottom bar below
     * the embedded SketcherWidget.
     */
    QSize minimumSizeHint() const override;
    QSize sizeHint() const override;

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

    /**
     * Allow this dialog's title bar to shrink so it's slightly shorter than the
     * atomistic/monomeric toggle buttons. That way, we can make the entire
     * dialog the same height as a standard SketcherWidget.
     */
    void configureWasmTitleBar();

    /**
     * Overridden Qt methods to update the bottom bar location when the dialog
     * is resized or laid out.
     */
    void resizeEvent(QResizeEvent* event) override;
    bool event(QEvent* event) override;

    std::unique_ptr<Ui::CustomMonomerDialog> ui;
    rdkit_extensions::ChainType m_chain_type;
    std::vector<int> m_required_attachment_points;

  private:
    /**
     * Calculate the dialog's size hint with the bottom bar either below the
     * SketcherWidget or in the SketcherWidget, but don't actually change the
     * footer's placement.
     */
    QSize dialogSizeHint(bool minimum, bool footer_below_view) const;

    /**
     * Place the bottom bar below the embedded SketcherWidget if there's enough
     * room, and move the bottom bar into the SketcherWidget if there's not.
     */
    void updateButtonBarPlacement();

    bool m_layout_ready = false;
    bool m_updating_layout = false;
    bool m_footer_below_view = false;
};

} // namespace sketcher
} // namespace schrodinger
