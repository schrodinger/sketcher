#pragma once

#include <memory>

#include <QPointer>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/monomer_database_filter.h"
#include "schrodinger/sketcher/widget/sketcher_view_with_wasm_outline_fix.h"

namespace Ui
{
class MonomerDatabaseFilterPopup;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * The Monomer Database popup's nested popup for editing its filter.
 *
 * The popup holds a copy of the filter. Every user change updates that copy
 * and immediately emits filterChanged(). The combo boxes never hold state of
 * their own: after every change, they're rewritten from the filter's
 * effective values, so an inapplicable (disabled) combo shows the most
 * permissive value, and a re-enabled one shows the value that the filter
 * remembered for it.
 */
class SKETCHER_API MonomerDatabaseFilterPopup
    : public SketcherViewWithWasmOutlineFix
{
    Q_OBJECT
  public:
    /**
     * @param anchor The widget that the popup is shown next to, i.e. the
     * Monomer Database popup's filter button
     * @param parent The parent widget
     */
    MonomerDatabaseFilterPopup(QWidget* anchor, QWidget* parent = nullptr);
    ~MonomerDatabaseFilterPopup();

    /**
     * Show the popup next to its anchor, initialized from the given filter
     */
    void showPopup(const MonomerDatabaseFilter& filter);

    /**
     * Move the open popup next to its anchor, e.g. after the anchor moved or
     * the available area changed. This does nothing if the popup is hidden.
     */
    void reposition();

    /**
     * @return The filter as edited in this popup
     */
    const MonomerDatabaseFilter& getFilter() const;

  signals:
    /**
     * Emitted whenever the user changes the filter
     */
    void filterChanged(const MonomerDatabaseFilter& filter);

  protected:
    /**
     * Overridden QWidget method: close any open combo box dropdown
     */
    void hideEvent(QHideEvent* event) override;

  private:
    /**
     * Respond to the user checking or unchecking a polymer type
     */
    void onPolymerTypesChanged();

    /**
     * Respond to the user choosing a monomer type
     */
    void onMonomerTypeChanged();

    /**
     * Respond to the user choosing a natural analog
     */
    void onNaturalAnalogChanged();

    /**
     * Reset the filter to the most permissive filter
     */
    void onClearAllClicked();

    /**
     * Update the widgets from the filter and emit filterChanged()
     */
    void onFilterEdited();

    /**
     * Write the filter's effective values, applicability and filtering state
     * to the widgets without emitting any signals
     */
    void updateWidgets();

    /**
     * Rebuild the natural analog combo box's items for the current polymer
     * types and select the filter's effective natural analog
     */
    void updateNaturalAnalogCombo();

    std::unique_ptr<Ui::MonomerDatabaseFilterPopup> m_ui;
    MonomerDatabaseFilter m_filter;
    // The Monomer Database popup's filter button
    QPointer<QWidget> m_anchor;
};

} // namespace sketcher
} // namespace schrodinger
