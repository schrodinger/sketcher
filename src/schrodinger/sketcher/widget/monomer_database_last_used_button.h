#pragma once

#include <QString>
#include <QToolButton>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h" // ChainType
#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A side bar button that shows the symbol of a page's last-used database
 * monomer. A long symbol is truncated and faded out the same way as the
 * monomer graphics items' labels.
 *
 * The checked state is controlled by the owning tool widget rather than by
 * clicks. Without a monomer, the button is blank and disabled, so that it
 * looks and acts like empty space.
 */
class SKETCHER_API MonomerDatabaseLastUsedButton : public QToolButton
{
    Q_OBJECT
  public:
    MonomerDatabaseLastUsedButton(QWidget* parent = nullptr);

    /**
     * Show the given monomer, or a blank button if its symbol is empty. Its
     * name is looked up in the monomer database for the tool tip. Nothing is
     * looked up if the monomer is the one that's already shown.
     */
    void setMonomer(const rdkit_extensions::MonomerID& monomer);

    /**
     * @return The monomer that the button shows, which has an empty symbol if
     * the button is blank
     */
    const rdkit_extensions::MonomerID& getMonomer() const;

    /**
     * Look up the shown monomer's name in the monomer database again, e.g.
     * after the database changes
     */
    void updateMonomerInfo();

    /**
     * @return Whether the symbol is too long to be shown in full, in which case
     * it's truncated and faded out
     */
    bool isTextFadedOut() const;

  protected:
    /**
     * Overridden QToolButton method: paint a truncated symbol with a fade-out
     */
    void paintEvent(QPaintEvent* event) override;

    /**
     * Overridden QToolButton method: preserve the model-controlled checked
     * state when the button is clicked
     */
    void nextCheckState() override;

  private:
    rdkit_extensions::MonomerID m_monomer{"",
                                          rdkit_extensions::ChainType::PEPTIDE};
    QString m_elided_text;
    QString m_fade_out_text;
};

} // namespace sketcher
} // namespace schrodinger
