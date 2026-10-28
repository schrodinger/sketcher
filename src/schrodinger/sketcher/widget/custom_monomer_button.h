#pragma once

#include <string>

#include <QToolButton>

#include "schrodinger/sketcher/definitions.h"

class QAction;
class QMenu;

namespace schrodinger
{
namespace rdkit_extensions
{
enum class ChainType;
}

namespace sketcher
{

/**
 * Button offering custom monomer drawing dialogs through an asynchronous menu.
 * Its checked state is controlled by the owning tool widget, rather than
 * clicks.
 */
class SKETCHER_API CustomMonomerButton : public QToolButton
{
    Q_OBJECT

  public:
    CustomMonomerButton(QWidget* parent = nullptr);

    /**
     * Set the first menu action to draw a peptide residue or nucleic acid base.
     */
    void setChainType(rdkit_extensions::ChainType chain_type);

  signals:
    /**
     * Forward the custom monomer accepted in a drawing dialog to the tool
     * widget.
     */
    void customMonomerAccepted(const std::string& smiles,
                               rdkit_extensions::ChainType monomer_type);

  protected:
    /**
     * Overriden QToolButton method: Paint centered contents with a standard
     * menu arrow on hover or check.
     */
    void paintEvent(QPaintEvent* event) override;

    /**
     * Overriden QToolButton method: Preserve the model-controlled checked state
     * when opening the menu.
     */
    void nextCheckState() override;

  private:
    QMenu* m_menu = nullptr;
    QAction* m_monomer_action = nullptr;
    rdkit_extensions::ChainType m_chain_type;

    /**
     * Open a drawing dialog and forward its accepted monomer.
     */
    void sketchCustomMonomer(rdkit_extensions::ChainType chain_type);
};

} // namespace sketcher
} // namespace schrodinger
