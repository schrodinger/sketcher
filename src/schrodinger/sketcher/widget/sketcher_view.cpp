#include "schrodinger/sketcher/widget/sketcher_view.h"

#include <QAction>
#include <QVariant>

#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

SetterPacket::SetterPacket(ModelKey key, std::function<void(int)> setter) :
    key(key),
    setter(setter)
{
}

SignalPacket::SignalPacket(ModelKey key, QAction* action) :
    key(key),
    action(action)
{
}

SketcherView::SketcherView(QWidget* parent) : QWidget(parent)
{
}

void SketcherView::setModel(SketcherModel* model)
{
    generatePackets();
    connectLocalSlots();
    if (m_sketcher_model == nullptr) {
        m_sketcher_model = model;
    } else {
        throw std::runtime_error(
            "The model has already been set on this view.");
    }
    connectToModel();
    updateWidgetsEnabled();
    updateCheckState();
}

SketcherModel* SketcherView::getModel() const
{
    return m_sketcher_model;
}

void SketcherView::updateWidgetsEnabled()
{
}

void SketcherView::updateCheckState()
{
}

std::vector<SetterPacket> SketcherView::getSetterPackets() const
{
    return m_setter_packets;
}

void SketcherView::onModelValuesChanged(
    const std::unordered_set<ModelKey>& keys)
{
    // In selection-only mode (e.g. in 2D Overlay) toolbars
    // are hidden; don't update widgets as there is a performance hit
    if (!parent() || isVisibleTo(window())) {
        updateWidgetsEnabled();
    }
    for (auto key : keys) {
        setValue(key, getModel()->getValue(key));
    }
    updateCheckState();
}

void SketcherView::onModelValuePinged(ModelKey key, QVariant value)
{
}

void SketcherView::setValue(ModelKey key, QVariant value)
{
    for (auto& setter_packet : getSetterPackets()) {
        if (key == setter_packet.key) {
            setter_packet.setter(value.toInt());
            return;
        }
    }
}

void SketcherView::generatePackets()
{
}

void SketcherView::connectLocalSlots()
{
    for (auto& signal_packet : m_signal_packets) {
        connect(signal_packet.action, &QAction::triggered, this,
                std::bind(&SketcherView::viewStateChanged, this,
                          signal_packet.key, std::placeholders::_1));
    }
}

void SketcherView::connectToModel()
{
    auto model = getModel();
    for (auto& setter_packet : getSetterPackets()) {
        QVariant value = model->getValue(setter_packet.key);
        setValue(setter_packet.key, value);
    }
    connect(this, &SketcherView::viewStateChanged, model,
            [model](auto key, auto value) { model->setValue(key, value); });
    connect(model, &SketcherModel::valuesChanged, this,
            &SketcherView::onModelValuesChanged);
    connect(model, &SketcherModel::valuePinged, this,
            &SketcherView::onModelValuePinged);
    connect(model, &SketcherModel::selectionChanged, this,
            &SketcherView::updateWidgetsEnabled);
    connect(model, &SketcherModel::selectionChanged, this,
            &SketcherView::updateCheckState);
    connect(model, &SketcherModel::interactiveItemsChanged, this,
            &SketcherView::updateWidgetsEnabled);
}

void SketcherView::disconnectUpdateWidgetsEnabled()
{
    auto model = getModel();
    disconnect(model, &SketcherModel::selectionChanged, this,
               &SketcherView::updateWidgetsEnabled);
    disconnect(model, &SketcherModel::interactiveItemsChanged, this,
               &SketcherView::updateWidgetsEnabled);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/sketcher_view.moc"
