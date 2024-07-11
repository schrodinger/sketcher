/* -------------------------------------------------------------------------
 * Tests class schrodinger::sketcher::EditAtomPropertiesDialog
 --------------------------------------------------------------------------- */

#pragma once

#include <memory>

#include <QValidator>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/modal_dialog.h"
#include "schrodinger/sketcher/rdkit/atom_properties.h"

class QLineEdit;
class QString;
class QTimer;
class QWidget;

namespace RDKit
{
class Atom;
}

namespace Ui
{
class EditAtomPropertiesDialog;
} // namespace Ui

namespace schrodinger
{
namespace sketcher
{

class MolModel;
class PeriodicTableWidget;

/**
 * Dialog used to edit all atom/query properties, including the ability to
 * toggle between setting as a fixed element atom or an atom query
 */
class SKETCHER_API EditAtomPropertiesDialog : public ModalDialog
{
  public:
    /**
     * @param atom The atom to be modified by the dialog
     * @param mol_model The MolModel containing the atom
     * @param parent The Qt parent widget
     */
    EditAtomPropertiesDialog(const RDKit::Atom* atom, MolModel* const mol_model,
                             QWidget* parent = nullptr);
    ~EditAtomPropertiesDialog();

    /**
     * Switch the dialog to the query page and the allowed list query type
     */
    void switchToAllowedList();

  protected:
    std::unique_ptr<Ui::EditAtomPropertiesDialog> ui;
    const RDKit::Atom* m_atom = nullptr;
    MolModel* m_mol_model = nullptr;
    QTimer* m_update_buttons_enabled_timer = nullptr;
    PeriodicTableWidget* m_atom_periodic_table_wdg = nullptr;
    PeriodicTableWidget* m_query_periodic_table_wdg = nullptr;

    // overridden QDialog method
    void accept() override;

    /**
     * Reset the dialog settings back to the input atom's values
     */
    void reset();

    /**
     * Replace the dark gray icon used for the clear button for all line edits
     * in the dialog with a light gray version
     */
    void changeAllLineEditClearButtonIcons();

    /**
     * Get the atom properties that are currently set in the dialog
     * @throw InvalidAtomPropertyError if the user has entered an invalid
     * element
     */
    std::shared_ptr<AbstractAtomProperties> getDialogSettings() const;

    /**
     * Get the SMARTS query entered by the user
     * @throw InvalidAtomPropertyError if the SMARTS query cannot be parsed
     */
    std::string getSmartsQuery() const;

    /**
     * Load the specified atom or atom query properties into the dialog
     */
    void loadProperties(const std::shared_ptr<AbstractAtomProperties> props);

    /**
     * Load the specified atom properties into the atom page of the dialog
     */
    void loadAtomProperties(const AtomProperties& props);

    /**
     * Load the specified atom query properties into the query page of the
     * dialog
     */
    void loadQueryProperties(const AtomQueryProperties& props);

    /**
     * Clear all of the query type entry fields (e.g. the specific element or
     * allowed list line edits) or reset them to their default values (e.g. the
     * wildcard combo box)
     */
    void clearQueryTypeFields();

    /**
     * Transfer any applicable properties between the atom page and the query
     * page whenever the user toggles between the two pages. E.g., copying the
     * atom element to the query element.
     * @param page_id The stacked widget ID of the page that was just switched
     * to
     */
    void onAtomOrQueryToggled(const int page_id);

    /**
     * Copy the contents of a single element line edit to the element list line
     * edit. This method will be a no-op if either of the following are true:
     *   - The single element line edit does not contain a valid atomic symbol
     *   - The element list line edit already starts with the contents of the
     *     single element line edit. If this happens, the single element line
     *     edit may have been populated by truncating the contents of the
     *     element list line edit. We don't want to truncate the element list
     *     line edit just because the user switched from query to atom and back,
     *     so we refrain from overwriting the element line contents in this
     *     scenario.
     * @param element_le The single element line edit to copy from
     */
    void
    transferSpecificElementToElementList(const QLineEdit* const element_le);

    /**
     * Copy the contents of the element list line edit to a single element line
     * edit. Only the first element of the element list will be transferred.
     * This method will be a no-op if the element list is empty, or if the first
     * item on the list is not a valid element.
     * @param element_le The single element line edit to copy to
     */
    void transferElementListToSpecificElement(QLineEdit* const element_le);

    /**
     * Update whether the OK and Reset buttons are enabled or disabled based on
     * the properties currently specified in the dialog. Both buttons will be
     * disabled if the current dialog properties are identical to the input
     * atom's properties. The OK button will also be disabled if the dialog
     * properties are invalid (e.g. an invalid atomic symbol).
     */
    void updateButtonsEnabled();

    /**
     * Update the dialog in response to a new selection in the query type combo
     * box.
     * @param combo_index The newly selected combo box index
     */
    void onQueryTypeComboBoxChanged(const int combo_index);

    /**
     * Update the atom (non-query) element line edit in response to a selection
     * in the periodic table widget
     * @param element The selected element
     */
    void onAtomPeriodicTableElementSelected(Element element);

    /**
     * Update the query element line edit in response to a selection in the
     * periodic table widget
     * @param element The selected element
     */
    void onQueryPeriodicTableElementSelected(Element element);
};

/**
 * A validator for line edits used to enter atomic symbols. The validator only
 * allows letters to be entered (as well as space and comma when allow_list is
 * true), and only reports the input as acceptable if it specifies real
 * elements. (E.g., "C" is acceptable input, "Qz" is not.)
 */
class ElementValidator : public QValidator
{
  public:
    /**
     * @param allow_list Whether the validator accepts a list of elements, or
     * only a single element. If true, commas and spaces will be allowed.
     * @param parent the Qt parent object
     */
    ElementValidator(const bool allow_list, QObject* parent = nullptr);
    QValidator::State validate(QString& input, int& pos) const override;

  protected:
    bool m_allow_list;
};

/**
 * An exception thrown when the user enters an invalid value for a property
 * (e.g. an atomic symbol that doesn't correspond to any element)
 */
class InvalidAtomPropertyError : public std::runtime_error
{
  public:
    using std::runtime_error::runtime_error;
};

} // namespace sketcher
} // namespace schrodinger
