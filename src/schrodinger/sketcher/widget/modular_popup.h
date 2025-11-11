#pragma once
#include <QToolButton>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view_with_wasm_outline_fix.h"

class QButtonGroup;
class QIcon;

namespace schrodinger
{
namespace sketcher
{

struct SKETCHER_API ButtonPacket {
    explicit ButtonPacket(QToolButton* button, int enum_int);
    QToolButton* button;
    int enum_int;
};

/**
 * A popup widget meant to hold a single button group.
 *
 * This class is meant to be used in conjunction with `ModularToolButton`.
 * Checking any of the buttons assigned to a concrete subclass of this widget is
 * meant to indicate a selection of a new appearance/functionality to the parent
 * `ModularToolButton` instance.
 */
class SKETCHER_API ModularPopup : public SketcherViewWithWasmOutlineFix
{
    Q_OBJECT

  public:
    ModularPopup(QWidget* parent = nullptr);

    /**
     * @param enum_int An enum item integer associated with a button in this
     * widget
     * @return The icon assigned to the button corresponding to `enum_int`. If
     * there is no such button, return a blank icon
     */
    QIcon getIcon(int enum_int) const;

    /**
     * @param enum_int An enum item integer associated with a button in this
     * widget
     * @return The text assigned to the button corresponding to `enum_int`. If
     * there is no such button, return the empty string
     */
    QString getText(int enum_int) const;

    /**
     * @param enum_int An enum item integer associated with a button in this
     * widget
     * @return The text assigned to the button corresponding to `enum_int`. If
     * there is no such button, return the empty string
     */
    QString getToolTip(int enum_int) const;

    /**
     * Return internal state of this widget to be used for testing.
     *
     * @return The button packets from this widget that are used to associate
     * buttons with enum item integers
     */
    std::vector<ButtonPacket> getButtonPackets() const;

    /**
     * @return IDs for all buttons
     */
    std::unordered_set<int> getButtonIDs() const;

  signals:
    void selectionChanged(int enum_int);

  protected:
    /**
     * Assign a button group instance for this widget.
     *
     * This method must be called exactly once by concrete subclasses to finish
     * initializing this widget.
     *
     * @param group The button group that contains all relevant buttons for this
     * widget
     */
    void setButtonGroup(QButtonGroup* group);

    /**
     * Generate and store button packets for this widget.
     *
     * This method should be implemented in concrete subclasses to populate
     * `m_button_packets` that associate enum item integers with the buttons on
     * this widget. The enum item integers indicate the button's intended
     * functionality.
     */
    virtual void generateButtonPackets() = 0;

    /**
     * Assign the specified tool button to be checked based on the model state.
     */
    void updateCheckState() override;

    /**
     * Get button ID to check.
     *
     * Meant to be reimplemented in concrete subclasses.
     *
     * @return The button ID that should be checked, based on the state of the
     * model. If no button should be checked, return -1.
     */
    virtual int getButtonIDToCheck() = 0;

    std::vector<ButtonPacket> m_button_packets;
    QButtonGroup* m_group = nullptr;

  protected slots:
    /**
     * Respond to the user clicking a button by emitting a signal and closing.
     *
     * @param button_id The group ID of the clicked button
     */
    void onButtonClicked(int button_id);
};

} // namespace sketcher
} // namespace schrodinger
