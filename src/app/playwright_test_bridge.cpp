/* -------------------------------------------------------------------------
 * Playwright-only test bridge for the WebAssembly Sketcher application.
 * ---------------------------------------------------------------------------
 */

#include "playwright_test_bridge.h"

#include <memory>
#include <stdexcept>
#include <unordered_map>

#include <QAbstractButton>
#include <QAction>
#include <QApplication>
#include <QClipboard>
#include <QComboBox>
#include <QEvent>
#include <QGraphicsView>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLineEdit>
#include <QMetaObject>
#include <QMenu>
#include <QPlainTextEdit>
#include <QPointer>
#include <QSpinBox>
#include <QTimer>
#include <QDoubleSpinBox>
#include <QWidget>

#ifdef __EMSCRIPTEN__
#include <emscripten/emscripten.h>
#endif

#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>

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

QJsonObject rect_to_object(const QPoint& top_left, const QRect& rect)
{
    QJsonObject result;
    result["x"] = top_left.x();
    result["y"] = top_left.y();
    result["width"] = rect.width();
    result["height"] = rect.height();
    return result;
}

std::string rect_to_json(const QPoint& top_left, const QRect& rect)
{
    return QJsonDocument(rect_to_object(top_left, rect))
        .toJson(QJsonDocument::Compact)
        .toStdString();
}

class MenuGeometryCache final : public QObject
{
  public:
    explicit MenuGeometryCache(SketcherWidget& sketcher) : m_sketcher(sketcher)
    {
        qApp->installEventFilter(this);
    }

    bool eventFilter(QObject* watched, QEvent* event) override
    {
        if (event->type() == QEvent::Show) {
            if (auto* widget = qobject_cast<QWidget*>(watched)) {
                cache_widget_geometry(*widget);
                const QPointer<QWidget> visible_widget(widget);
                QTimer::singleShot(0, this, [this, visible_widget] {
                    if (visible_widget) {
                        cache_widget_geometry(*visible_widget);
                    }
                });
            }
            if (auto* menu = qobject_cast<QMenu*>(watched)) {
                cache_menu_geometry(*menu);
            }
        }
        return false;
    }

    std::string action_rect(const std::string& object_name_or_text) const
    {
        const auto query = QString::fromStdString(object_name_or_text);
        // A menu may first be shown at a provisional location while Qt lays
        // out a cascading submenu.  The event-time cache is useful for
        // publishing diagnostics, but must not be trusted for an interaction:
        // find the action in the menu currently visible on screen instead.
        for (auto* widget : QApplication::allWidgets()) {
            const auto* menu = qobject_cast<QMenu*>(widget);
            if (!menu || !menu->isVisible()) {
                continue;
            }
            for (auto* action : menu->actions()) {
                if (action->objectName() != query && action->text() != query) {
                    continue;
                }
                const QRect action_rect = menu->actionGeometry(action);
                if (!action_rect.isEmpty()) {
                    return rect_to_json(
                        menu->mapTo(&m_sketcher, action_rect.topLeft()),
                        action_rect);
                }
            }
        }
        const auto result = m_action_rects.find(object_name_or_text);
        return result == m_action_rects.end() ? "{}" : result->second;
    }

  private:
    void cache_menu_geometry(const QMenu& menu)
    {
        QJsonObject browser_rects;
        for (auto* action : menu.actions()) {
            const QRect action_rect = menu.actionGeometry(action);
            if (action_rect.isEmpty()) {
                continue;
            }
            const auto result = rect_to_json(
                menu.mapTo(&m_sketcher, action_rect.topLeft()), action_rect);
            const auto browser_result = rect_to_object(
                menu.mapTo(&m_sketcher, action_rect.topLeft()), action_rect);
            if (!action->objectName().isEmpty()) {
                m_action_rects[action->objectName().toStdString()] = result;
                browser_rects[action->objectName()] = browser_result;
            }
            m_action_rects[action->text().toStdString()] = result;
            browser_rects[action->text()] = browser_result;
        }
        publish_browser_rects(browser_rects);
    }

    void publish_browser_rects(const QJsonObject& rects) const
    {
#ifdef __EMSCRIPTEN__
        const auto json = QJsonDocument(rects).toJson(QJsonDocument::Compact);
        EM_ASM(
            {
                window.__sketcherPlaywrightMenuRects =
                    JSON.parse(UTF8ToString($0));
            },
            json.constData());
#endif
    }

    void cache_widget_geometry(const QWidget& widget)
    {
        cache_widget(widget);
        for (auto* child : widget.findChildren<QWidget*>()) {
            cache_widget(*child);
        }
        publish_browser_widget_rects();
    }

    void cache_widget(const QWidget& widget)
    {
        if (!widget.isVisible() || widget.objectName().isEmpty()) {
            return;
        }
        m_widget_rects[widget.objectName()] = rect_to_object(
            widget.mapTo(&m_sketcher, QPoint(0, 0)), widget.rect());
    }

    void publish_browser_widget_rects() const
    {
#ifdef __EMSCRIPTEN__
        const auto json =
            QJsonDocument(m_widget_rects).toJson(QJsonDocument::Compact);
        EM_ASM(
            {
                window.__sketcherPlaywrightWidgetRects =
                    JSON.parse(UTF8ToString($0));
            },
            json.constData());
#endif
    }

    SketcherWidget& m_sketcher;
    std::unordered_map<std::string, std::string> m_action_rects;
    QJsonObject m_widget_rects;
};

MenuGeometryCache& menu_geometry_cache(SketcherWidget& sketcher)
{
    static auto cache = std::make_unique<MenuGeometryCache>(sketcher);
    return *cache;
}

QGraphicsView& require_view(SketcherWidget& sketcher)
{
    return dynamic_cast<QGraphicsView&>(
        require_visible_widget(sketcher, "view"));
}

std::string scene_item_rect(SketcherWidget& sketcher, const QGraphicsItem& item)
{
    auto& view = require_view(sketcher);
    const auto viewport_rect =
        view.mapFromScene(item.sceneBoundingRect()).boundingRect();
    const auto top_left =
        view.viewport()->mapTo(&sketcher, viewport_rect.topLeft());
    return rect_to_json(top_left, viewport_rect);
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
    // Install the menu observer before Playwright opens the first menu.
    menu_geometry_cache(sketcher);
    auto* widget =
        find_visible_widget(sketcher, QString::fromStdString(object_name));
    if (!widget) {
        return "{}";
    }
    const QPoint top_left = widget->mapTo(&sketcher, QPoint(0, 0));
    return rect_to_json(top_left, widget->rect());
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

std::string get_menu_action_rect(SketcherWidget& sketcher,
                                 const std::string& object_name_or_text)
{
    return menu_geometry_cache(sketcher).action_rect(object_name_or_text);
}

std::string get_atom_rect(SketcherWidget& sketcher, const int atom_index)
{
    for (auto* item : require_view(sketcher).scene()->items()) {
        auto* atom_item = qgraphicsitem_cast<AtomItem*>(item);
        if (atom_item && static_cast<int>(atom_item->getAtom()->getIdx()) + 1 ==
                             atom_index) {
            return scene_item_rect(sketcher, *atom_item);
        }
    }
    return "{}";
}

std::string get_bond_rect(SketcherWidget& sketcher, const int bond_index)
{
    for (auto* item : require_view(sketcher).scene()->items()) {
        auto* bond_item = qgraphicsitem_cast<BondItem*>(item);
        if (bond_item && static_cast<int>(bond_item->getBond()->getIdx()) + 1 ==
                             bond_index) {
            return scene_item_rect(sketcher, *bond_item);
        }
    }
    return "{}";
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
