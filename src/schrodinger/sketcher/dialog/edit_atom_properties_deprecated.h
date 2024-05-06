/* -------------------------------------------------------------------------
 * Tests class schrodinger::sketcher::EditAtomPropertiesDialog
 --------------------------------------------------------------------------- */

#pragma once

#include <memory>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/dialog/modal_dialog.h"
#include "schrodinger/sketcher/rdkit/atom_properties.h"

class QComboBox;
class QLineEdit;
class QRadioButton;
class QSpinBox;
class sketcherAtom;
enum class EnhancedStereoType;

namespace Ui
{
class EditAtomPropertiesDialogDeprecated;
class CommonAtomPropertiesWidget;
} // namespace Ui

namespace schrodinger
{
namespace sketcher
{

enum class AtomQuery;
enum class Element;
class PeriodicTableWidget;
class SketcherModel;

/**
 * Dialog used to edit all atom/query properties, including the ability to
 * toggle between setting as a fixed element atom or an atom query
 */
class SKETCHER_API EditAtomPropertiesDialogDeprecated : public ModalDialog
{
    Q_OBJECT
  public:
    EditAtomPropertiesDialogDeprecated(SketcherModel* model, sketcherAtom& atom,
                                       QWidget* parent = nullptr,
                                       Qt::WindowFlags f = Qt::WindowFlags());
    ~EditAtomPropertiesDialogDeprecated();

    /**
     * Updates the dialog to assign "Allowed List" as the option for "Element".
     */
    void setToAllowedList();

  protected slots:
    void accept() override;

    /**
     * Reset the widget to the current atom values on the atom
     */
    void reset();

    /**
     * Organizational methods to reset different components of the dialog. Only
     * meant to be called from the `reset()` method.
     */
    void resetGeneralQuerySubwidgets();
    void resetAdvancedQuerySubwidgets();

    /**
     * Writes values from the atom properties page to the member sketcher atom
     */
    void writeAtomInfo();

    /**
     * Writes values from the query properties page to the member sketcher atom
     */
    void writeQueryInfo();

    /**
     * Respond when the "query type" combo box changes.
     *
     * Specifically, update properties of the periodic table button/widget.
     */
    void onQueryTypeComboValueChanged();

    /**
     * Respond when the user clicks a button in the periodic table widget.
     *
     * Either add the selected element to the existing query type text or
     * replace the text, depending on the mode.
     *
     * @param element The clicked element
     */
    void onPeriodicTableElementSelected(Element element);

    /**
     * Enable or disable OK button based on state of dialog.
     */
    void updateOKButtonEnabled();

  protected:
    std::unique_ptr<Ui::EditAtomPropertiesDialogDeprecated> ui;
    sketcherAtom& m_atom;
    QueryType getQueryTypeComboValue() const;
    AtomQuery getWildcardComboValue() const;
    void setQueryTypeComboValue(QueryType etype);
    PeriodicTableWidget* m_periodic_table_wdg = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
};

/**
 * Widget editing common properties between atoms and queries, including
 * isotope, charge, number of unpaired electrons, and enhances stereo labels.
 * Used specifically in the context of EditAtomPropertiesDialog
 */
class SKETCHER_API CommonAtomPropertiesWidget : public QWidget
{
    Q_OBJECT
  public:
    CommonAtomPropertiesWidget(QWidget* parent = nullptr);
    ~CommonAtomPropertiesWidget();

    /**
     * Populates the widget with properties read from the given atom
     * @param atom sketcher atom from which to read
     */
    void readAtomInfo(const sketcherAtom& atom);

    /**
     * Extracts values from the widget and sets them as atom properties
     * @param atom sketcher atom to which to write
     */
    void writeAtomInfo(sketcherAtom& atom) const;

    /**
     * Update enabled elements when the main dialog query combo changes.
     * @param query_type current "query type" combo box value
     */
    void onQueryTypeComboValueChanged(QueryType query_type);

    /**
     * @param enable Whether to enable the enhanced stereo combo box
     */
    void setEnhancedStereoTypeComboEnabled(bool enable);
    std::unique_ptr<Ui::CommonAtomPropertiesWidget> ui;

  protected slots:
    /**
     * Updates the enhanced stereo so that Off/ABS disable/clear the group ID,
     * whereas AND/OR enable and populate it
     */
    void updateEnhancedStereoSpinBox();

  protected:
    friend EditAtomPropertiesDialogDeprecated::
        EditAtomPropertiesDialogDeprecated(SketcherModel* model,
                                           sketcherAtom& atom, QWidget* parent,
                                           Qt::WindowFlags f);

    void setEnhancedStereoTypeComboValue(
        const ::EnhancedStereoType& enhanced_stereo_type);
};

} // namespace sketcher
} // namespace schrodinger
