#pragma once
#include <functional>
#include <unordered_set>

#include <QWidget>

#include "schrodinger/sketcher/definitions.h"

class QVariant;
class QAction;
namespace schrodinger
{
namespace sketcher
{

class SketcherModel;
enum class ModelKey;

/**
 * Data class to hold information about view setters.
 */
struct SetterPacket {
    explicit SetterPacket(ModelKey key, std::function<void(int)> setter);
    ModelKey key;
    std::function<void(int)> setter;
};

/**
 * Data class to hold information about view signals.
 */
struct SignalPacket {
    SignalPacket(ModelKey key, QAction* action);
    ModelKey key;
    QAction* action;
};

class SKETCHER_API SketcherView : public QWidget
{
    Q_OBJECT

  public:
    SketcherView(QWidget* parent = nullptr);

    ~SketcherView() = default;

    /**
     * Assign a model to this view.
     *
     * On assignment, the view should be synchronized to the state of the model
     * and all relevant signal/slot connections should be made. Note that the
     * model on a view should be set exactly once.
     *
     * @throw std::runtime_error If a model is set on this view more than once.
     */
    virtual void setModel(SketcherModel* model);

    /**
     * @return The model assigned to this view.
     */
    SketcherModel* getModel() const;

    /**
     * Get the setter packets used to store model signal/view slot information.
     *
     * This method is public only so that it can be used for unit testing.
     */
    std::vector<SetterPacket> getSetterPackets() const;

    /**
     * Assign the value of some model-mapped state to this view.
     *
     * @param key The model key associated with the state
     * @param value The value stored on the model for this state
     */
    void setValue(ModelKey key, QVariant value);

    /**
     * Disconnect signals from changes to the model to trigger
     * calls to this view's updateWidgetsEnabled() slot.
     */
    void disconnectUpdateWidgetsEnabled();

  signals:
    void viewStateChanged(ModelKey key, QVariant value);

  protected:
    /**
     * Data necessary for updating the state of the model from checkable
     * actions on this view.
     *
     * These are used when synchronizing the model to this view.
     */
    std::vector<SignalPacket> m_signal_packets;

    /**
     * Data necessary for updating the state of this view from the model.
     *
     * These are used when synchronizing this view to the model.
     */
    std::vector<SetterPacket> m_setter_packets;

    /**
     * Generate a vector that identifies setters associated with each model key.
     *
     * This function is meant to be run a single time during instantiation in
     * order to assign relevant values to `m_setter_packets` and
     * `m_signal_packets`.
     */
    virtual void generatePackets();

    /**
     * Connect signals and slots owned by this widget (or its subwidgets).
     *
     * This function is meant to be called a single time after
     * `generatePackets()`, as it connects signals and slots defined in signal
     * packets. It does not need to be called by subclasses, as it will be
     * called in the default implementation of `setModel()`.
     */
    virtual void connectLocalSlots();

    /**
     * Connect model signals/slots and synchronize view to current model state.
     */
    virtual void connectToModel();

  protected slots:
    /**
     * Respond to changes in the model's value map by updating this view.
     *
     * @param keys The keys associated with values that have been changed in the
     * model
     */
    virtual void onModelValuesChanged(const std::unordered_set<ModelKey>& keys);

    /**
     * Respond to a value being pinged on the model.
     *
     * A value is pinged when it gets set, but can also be pinged directly
     * without attempting to change the value of the model.
     *
     * This slot doesn't do anything, but it should be overridden by subclasses.
     *
     * @param key The key of the pinged value
     * @param value The value that the key was pinged with
     */
    virtual void onModelValuePinged(ModelKey key, QVariant value);

    /**
     * Update which widgets are enabled depending on the state of the model.
     *
     * Should be overridden in subclasses that enable/disable certain widgets
     * when the model changes.
     */
    virtual void updateWidgetsEnabled();

    /**
     * Update tool button check states depending on the state of the model.
     *
     * Should be overridden in subclasses that dead with check states.
     */
    virtual void updateCheckState();

  private:
    /**
     * The sketcher model instance.
     *
     * Any state held by this view should be synchronized with the model.
     */
    SketcherModel* m_sketcher_model = nullptr;
};

} // namespace sketcher
} // namespace schrodinger
