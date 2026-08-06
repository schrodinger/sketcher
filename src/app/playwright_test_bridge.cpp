/* -------------------------------------------------------------------------
 * Playwright-only test bridge for the WebAssembly Sketcher application.
 * ---------------------------------------------------------------------------
 */

#include "playwright_test_bridge.h"

#include <stdexcept>

#include <QAbstractButton>
#include <QAction>
#include <QApplication>
#include <QClipboard>
#include <QComboBox>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLineEdit>
#include <QMetaObject>
#include <QPlainTextEdit>
#include <QSpinBox>
#include <QTimer>
#include <QDoubleSpinBox>
#include <QWidget>

#include "schrodinger/sketcher/sketcher_widget.h"

namespace schrodinger::sketcher::playwright_test_bridge
{
namespace
{
QWidget* find_visible_widget(SketcherWidget& sketcher, const QString& name)
{
    const auto candidates = sketcher.findChildren<QWidget*>(name);
    for (auto* widget : candidates) {
        if (widget->isVisible()) {
            return widget;
        }
    }
    return nullptr;
}

QWidget& require_visible_widget(SketcherWidget& sketcher,
                                const std::string& object_name)
{
    auto* widget =
        find_visible_widget(sketcher, QString::fromStdString(object_name));
    if (!widget) {
        throw std::runtime_error("playwright test bridge: visible widget not "
                                 "found with objectName '" +
                                 object_name + "'");
    }
    return *widget;
}

void set_text(QWidget& widget, const QString& text)
{
    if (auto* line_edit = qobject_cast<QLineEdit*>(&widget)) {
        line_edit->setText(text);
        return;
    }
    if (auto* text_edit = qobject_cast<QPlainTextEdit*>(&widget)) {
        text_edit->setPlainText(text);
        return;
    }
    if (auto* combo_box = qobject_cast<QComboBox*>(&widget)) {
        const int index = combo_box->findText(text);
        if (index < 0) {
            throw std::runtime_error("playwright test bridge: combo-box item "
                                     "not found: '" +
                                     text.toStdString() + "'");
        }
        combo_box->setCurrentIndex(index);
        return;
    }
    if (auto* spin_box = qobject_cast<QSpinBox*>(&widget)) {
        spin_box->setValue(text.toInt());
        return;
    }
    if (auto* spin_box = qobject_cast<QDoubleSpinBox*>(&widget)) {
        spin_box->setValue(text.toDouble());
        return;
    }
    throw std::runtime_error("playwright test bridge: widget does not support "
                             "text/value input");
}
} // namespace

void click_button(SketcherWidget& sketcher, const std::string& object_name)
{
    auto* button = sketcher.findChild<QAbstractButton*>(
        QString::fromStdString(object_name));
    if (!button) {
        throw std::runtime_error("playwright test bridge: button not found "
                                 "with objectName '" +
                                 object_name + "'");
    }
    button->click();
}

std::string get_widget_rect(SketcherWidget& sketcher,
                            const std::string& object_name)
{
    auto* widget =
        find_visible_widget(sketcher, QString::fromStdString(object_name));
    if (!widget) {
        return "{}";
    }
    const QPoint top_left = widget->mapTo(&sketcher, QPoint(0, 0));
    QJsonObject result;
    result["x"] = top_left.x();
    result["y"] = top_left.y();
    result["width"] = widget->width();
    result["height"] = widget->height();
    return QJsonDocument(result).toJson(QJsonDocument::Compact).toStdString();
}

std::string get_widget_state(SketcherWidget& sketcher,
                             const std::string& object_name)
{
    auto& widget = require_visible_widget(sketcher, object_name);
    QJsonObject result;
    result["objectName"] = widget.objectName();
    result["className"] = widget.metaObject()->className();
    result["visible"] = widget.isVisible();
    result["enabled"] = widget.isEnabled();

    if (auto* button = qobject_cast<QAbstractButton*>(&widget)) {
        result["checked"] = button->isChecked();
        result["text"] = button->text();
    } else if (auto* line_edit = qobject_cast<QLineEdit*>(&widget)) {
        result["text"] = line_edit->text();
    } else if (auto* text_edit = qobject_cast<QPlainTextEdit*>(&widget)) {
        result["text"] = text_edit->toPlainText();
    } else if (auto* combo_box = qobject_cast<QComboBox*>(&widget)) {
        result["text"] = combo_box->currentText();
    } else if (auto* spin_box = qobject_cast<QSpinBox*>(&widget)) {
        result["text"] = spin_box->text();
        result["value"] = spin_box->value();
    } else if (auto* spin_box = qobject_cast<QDoubleSpinBox*>(&widget)) {
        result["text"] = spin_box->text();
        result["value"] = spin_box->value();
    }
    return QJsonDocument(result).toJson(QJsonDocument::Compact).toStdString();
}

void set_widget_text(SketcherWidget& sketcher, const std::string& object_name,
                     const std::string& text)
{
    set_text(require_visible_widget(sketcher, object_name),
             QString::fromStdString(text));
}

void activate_menu_action(SketcherWidget& sketcher,
                          const std::string& object_name_or_text)
{
    const QString query = QString::fromStdString(object_name_or_text);
    const auto actions = sketcher.findChildren<QAction*>();
    for (auto* action : actions) {
        if (action->isVisible() &&
            (action->objectName() == query || action->text() == query)) {
            // An embind call runs while Qt is handling the browser event.
            // Triggering an action synchronously can re-enter Qt's event loop
            // and prevent the browser-side call from returning. Queue it for
            // the next Qt event-loop iteration instead.
            QTimer::singleShot(0, action, [action] { action->trigger(); });
            return;
        }
    }
    throw std::runtime_error("playwright test bridge: visible menu action not "
                             "found with objectName/text '" +
                             object_name_or_text + "'");
}

std::string clipboard_text()
{
    return QApplication::clipboard()->text().toStdString();
}

void set_clipboard_text(const std::string& text)
{
    QApplication::clipboard()->setText(QString::fromStdString(text));
}
} // namespace schrodinger::sketcher::playwright_test_bridge
