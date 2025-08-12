#pragma once
#include <memory>

#include <boost/bimap.hpp>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/abstract_draw_tool_widget.h"

class QAbstractButton;

namespace Ui
{
class SetAtomWidget;
}

namespace schrodinger
{
namespace sketcher
{

class AtomQueryPopup;
class PeriodicTableWidget;
enum class Element;

/**
 * Interface for user to select various atom tools, including queries.
 */
class SKETCHER_API SetAtomWidget : public AbstractDrawToolWidget
{
    Q_OBJECT

  public:
    SetAtomWidget(QWidget* parent = nullptr);
    ~SetAtomWidget();

    /**
     * Assign a model to this view.
     *
     * On assignment, the view should be synchronized to the state of the model
     * and all relevant signal/slot connections should be made. Note that the
     * model on a view should be set exactly once.
     *
     * @throw std::runtime_error If a model is set on this view more than once.
     */
    void setModel(SketcherModel* model) override;

  protected:
    std::unique_ptr<Ui::SetAtomWidget> ui;
    AtomQueryPopup* m_atom_query_wdg = nullptr;
    PeriodicTableWidget* m_periodic_table_wdg = nullptr;
    boost::bimap<QAbstractButton*, Element> m_button_element_bimap;
    void updateWidgetsEnabled() override;
    void updateCheckedButton() override;
    std::unordered_set<QAbstractButton*> getCheckableButtons() override;

    /**
     * @param button A button in this widget associated with an element (not a
     * query)
     * @return the element associated with the specified button
     */
    Element getElementForButton(QAbstractButton* button) const;

  protected slots:
    /**
     * Respond to a value being pinged on the model.
     *
     * Specifically, update the "last picked element" button by updating its
     * element value if an element is ever set on the model.
     */
    void onModelValuePinged(ModelKey key, QVariant value) override;

    void onAtomButtonClicked(QAbstractButton* button);
};

/**
 * Extension of set atom widget to be used as a pop-out in context menus
 */
class SKETCHER_API SetAtomMenuWidget : public SetAtomWidget
{
    Q_OBJECT

  public:
    SetAtomMenuWidget(QWidget* parent = nullptr);
    ~SetAtomMenuWidget();

  signals:
    void anyButtonClicked();
};

} // namespace sketcher
} // namespace schrodinger
