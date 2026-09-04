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
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QLabel>
#include <QLineEdit>
#include <QMouseEvent>
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
#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"
#include "schrodinger/sketcher/molviewer/abstract_bond_or_connector_item.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/widget/tool_button_with_popup.h"

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>

namespace schrodinger::sketcher::playwright_test_bridge
{
namespace
{
QWidget* find_visible_widget(SketcherWidget& sketcher, const QString& name)
{
    auto menu_ancestor = [](QWidget* widget) {
        for (auto* parent = widget ? widget->parentWidget() : nullptr; parent;
             parent = parent->parentWidget()) {
            if (qobject_cast<QMenu*>(parent)) {
                return true;
            }
        }
        return false;
    };
    // A modal dialog can legitimately share an objectName with the persistent
    // toolbar (for example, both have an export_btn).  A person can only
    // interact with the active dialog, so prefer its visible control before
    // considering the background Sketcher hierarchy.
    if (auto* modal = QApplication::activeModalWidget(); modal &&
        modal->isVisible()) {
        if (modal->objectName() == name) {
            return modal;
        }
        const auto modal_candidates = modal->findChildren<QWidget*>(name);
        for (auto* widget : modal_candidates) {
            if (widget->isVisible()) {
                return widget;
            }
        }
    }
    // A context-menu widget can legitimately share an objectName with the
    // persistent toolbar (for example periodic_table_btn). Prefer the visible
    // control contained in the open QMenu, so browser input follows the
    // source's context-menu path rather than silently using the toolbar.
    QWidget* fallback = nullptr;
    const auto candidates = sketcher.findChildren<QWidget*>(name);
    for (auto* widget : candidates) {
        if (!widget->isVisible()) {
            continue;
        }
        if (menu_ancestor(widget)) {
            return widget;
        }
        if (!fallback) {
            fallback = widget;
        }
    }
    // Modular tool popups are separate top-level Qt/WASM windows rather than
    // children of SketcherWidget.  Search those windows too, solely for their
    // visible named controls used by browser tests.
    for (auto* top_level : QApplication::topLevelWidgets()) {
        if (!top_level) {
            continue;
        }
        if (top_level->objectName() == name && top_level->isVisible()) {
            return top_level;
        }
        const auto popup_candidates = top_level->findChildren<QWidget*>(name);
        for (auto* widget : popup_candidates) {
            if (!widget->isVisible()) {
                continue;
            }
            if (menu_ancestor(widget)) {
                return widget;
            }
            if (!fallback) {
                fallback = widget;
            }
        }
    }
    return fallback;
}

QWidget* find_widget(SketcherWidget& sketcher, const QString& name)
{
    if (auto* widget = find_visible_widget(sketcher, name)) {
        return widget;
    }
    if (auto* widget = sketcher.findChild<QWidget*>(name)) {
        return widget;
    }
    for (auto* top_level : QApplication::topLevelWidgets()) {
        if (!top_level) {
            continue;
        }
        if (top_level->objectName() == name) {
            return top_level;
        }
        if (auto* widget = top_level->findChild<QWidget*>(name)) {
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

QWidget& require_widget(SketcherWidget& sketcher,
                        const std::string& object_name)
{
    auto* widget = find_widget(sketcher, QString::fromStdString(object_name));
    if (!widget) {
        throw std::runtime_error("playwright test bridge: widget not found "
                                 "with objectName '" +
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

QString normalise_menu_action_name(QString name)
{
    // QAction text may contain an ampersand to mark its mnemonic (for
    // example, "Modify &All").  The browser-facing wrappers intentionally
    // use the human-visible name, so compare the two forms consistently.
    name.remove('&');
    return name.simplified();
}

std::string rect_to_json(const QPoint& top_left, const QRect& rect)
{
    return QJsonDocument(rect_to_object(top_left, rect))
        .toJson(QJsonDocument::Compact)
        .toStdString();
}

QPoint map_popup_point_to_sketcher(const QWidget& popup,
                                   const SketcherWidget& sketcher,
                                   const QPoint& point)
{
    // A cascading QMenu is its own popup window.  Mapping it directly to the
    // Sketcher widget is fragile in Qt/WASM; convert through global
    // coordinates, as Qt does for independent top-level widgets.
    return sketcher.mapFromGlobal(popup.mapToGlobal(point));
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
                if (!qobject_cast<QMenu*>(widget)) {
                    cache_widget_geometry(*widget);
                    const QPointer<QWidget> visible_widget(widget);
                    QTimer::singleShot(0, this, [this, visible_widget] {
                        if (visible_widget) {
                            cache_widget_geometry(*visible_widget);
                        }
                    });
                }
            }
            // QMenu inspection is intentionally disabled pending a WASM-safe
            // geometry mechanism.  Accessing popup menu internals during its
            // show lifecycle aborts the current Qt/WASM runtime.
        }
        return false;
    }

    std::string action_rect(const std::string& object_name_or_text) const
    {
        const auto query = QString::fromStdString(object_name_or_text);
        const auto normalised_query = normalise_menu_action_name(query);
        const auto action_rect_in_menu = [&](QMenu* menu) -> std::string {
            if (!menu || !menu->isVisible()) {
                return {};
            }
            for (auto* action : menu->actions()) {
                if (action->objectName() != query && action->text() != query &&
                    normalise_menu_action_name(action->objectName()) !=
                        normalised_query &&
                    normalise_menu_action_name(action->text()) !=
                        normalised_query) {
                    continue;
                }
                const QRect action_rect = menu->actionGeometry(action);
                if (!action_rect.isEmpty()) {
                    return rect_to_json(
                        map_popup_point_to_sketcher(*menu, m_sketcher,
                                                    action_rect.topLeft()),
                        action_rect);
                }
            }
            return {};
        };
        // This is deliberately an on-demand lookup: inspecting QMenu while it
        // is being shown is unsafe in Qt/WASM, but once the popup is visible
        // its action geometry can be read for a real browser mouse click.
        // A nested context path can contain duplicate labels.  Qt identifies
        // the deepest currently open submenu as the active popup, which is
        // exactly the menu receiving the next human pointer event.
        if (auto* active_menu =
                qobject_cast<QMenu*>(QApplication::activePopupWidget())) {
            if (const auto result = action_rect_in_menu(active_menu);
                !result.empty()) {
                return result;
            }
        }
        for (auto* menu : m_sketcher.findChildren<QMenu*>()) {
            if (const auto result = action_rect_in_menu(menu);
                !result.empty()) {
                return result;
            }
        }
        const auto result = m_action_rects.find(object_name_or_text);
        if (result != m_action_rects.end()) {
            return result->second;
        }
        for (const auto& [name, rect] : m_action_rects) {
            if (normalise_menu_action_name(QString::fromStdString(name)) ==
                normalised_query) {
                return rect;
            }
        }
        return "{}";
    }

  private:
    void cache_menu_geometry(const QMenu& menu)
    {
        for (auto* action : menu.actions()) {
            const QRect action_rect = menu.actionGeometry(action);
            if (action_rect.isEmpty()) {
                continue;
            }
            const auto result =
                rect_to_json(map_popup_point_to_sketcher(menu, m_sketcher,
                                                         action_rect.topLeft()),
                             action_rect);
            if (!action->objectName().isEmpty()) {
                m_action_rects[action->objectName().toStdString()] = result;
            }
            m_action_rects[action->text().toStdString()] = result;
        }
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

std::string generic_item_rect(SketcherWidget& sketcher, const bool is_atom,
                              const int index)
{
    for (auto* item : require_view(sketcher).scene()->items()) {
        if (!item->isVisible()) {
            continue;
        }
        if (is_atom) {
            auto* atom_item = dynamic_cast<AbstractAtomOrMonomerItem*>(item);
            if (atom_item != nullptr &&
                static_cast<int>(atom_item->getAtom()->getIdx()) == index) {
                return scene_item_rect(sketcher, *atom_item);
            }
        } else {
            auto* bond_item = dynamic_cast<AbstractBondOrConnectorItem*>(item);
            if (bond_item != nullptr &&
                static_cast<int>(bond_item->getBond()->getIdx()) == index) {
                return scene_item_rect(sketcher, *bond_item);
            }
        }
    }
    return "{}";
}
} // namespace

std::string get_rect(SketcherWidget& sketcher, const std::string& selector)
{
    const auto separator = selector.find(':');
    if (separator == std::string::npos) {
        throw std::runtime_error("playwright test bridge: malformed selector '" +
                                 selector + "' (expected '<kind>:<value>')");
    }
    const auto kind = selector.substr(0, separator);
    const auto value = selector.substr(separator + 1);
    if (kind == "widget") {
        const auto rect = get_widget_rect(sketcher, value);
        if (rect == "{}") {
            return rect;
        }
        // Generic geometry consumers use enabled to avoid clicking a painted
        // control that Qt would reject. The older rectangle API predates that
        // field, so add it here without changing the legacy payload.
        auto object = QJsonDocument::fromJson(
                          QByteArray::fromStdString(rect))
                          .object();
        object["enabled"] =
            find_visible_widget(sketcher, QString::fromStdString(value))
                ->isEnabled();
        return QJsonDocument(object).toJson(QJsonDocument::Compact).toStdString();
    }
    if (kind == "state") {
        auto object = QJsonDocument::fromJson(
                          QByteArray::fromStdString(
                              get_widget_state(sketcher, value)))
                          .object();
        auto& widget = require_widget(sketcher, value);
        const QPoint top_left =
            sketcher.mapFromGlobal(widget.mapToGlobal(QPoint(0, 0)));
        object["x"] = top_left.x();
        object["y"] = top_left.y();
        object["width"] = widget.width();
        object["height"] = widget.height();
        return QJsonDocument(object).toJson(QJsonDocument::Compact).toStdString();
    }
    if (kind == "menu") {
        // MenuGeometryCache performs its lookup against the active visible
        // popup. Walking the full QAction hierarchy here just to add an
        // `enabled` field is unsafe while Qt/WASM is creating a cascading
        // submenu, and is not needed to obtain the browser click geometry.
        return get_menu_action_rect(sketcher, value);
    }
    if (kind == "atom" || kind == "bond") {
        size_t parsed = 0;
        int index = 0;
        try {
            index = std::stoi(value, &parsed);
        } catch (const std::exception&) {
            throw std::runtime_error("playwright test bridge: invalid " + kind +
                                     " index '" + value + "'");
        }
        if (parsed != value.size() || index < 0) {
            throw std::runtime_error("playwright test bridge: invalid " + kind +
                                     " index '" + value + "'");
        }
        // Generic selectors use the model's native zero-based atom/bond
        // indexes and cover monomer scene items as well as ordinary atoms.
        return generic_item_rect(sketcher, kind == "atom", index);
    }
    throw std::runtime_error("playwright test bridge: unrecognized selector kind '" +
                             kind + "'");
}

std::string get_popup_owner(SketcherWidget& sketcher,
                            const std::string& object_name)
{
    const auto wanted = QString::fromStdString(object_name);
    for (auto* button : sketcher.findChildren<ToolButtonWithPopup*>()) {
        auto* popup = button->getPopupWidget();
        if (popup != nullptr &&
            (popup->objectName() == wanted ||
             popup->findChild<QWidget*>(wanted) != nullptr)) {
            return button->objectName().toStdString();
        }
    }
    return {};
}

void activate_action(SketcherWidget& sketcher,
                     const std::string& object_name_or_text)
{
    const auto query = QString::fromStdString(object_name_or_text);
    const auto normalised_query = normalise_menu_action_name(query);
    for (auto* action : sketcher.findChildren<QAction*>()) {
        if (action->objectName() != query &&
            normalise_menu_action_name(action->text()) != normalised_query) {
            continue;
        }
        if (!action->isEnabled()) {
            throw std::runtime_error("playwright test bridge: action '" +
                                     object_name_or_text + "' is disabled");
        }
        action->trigger();
        return;
    }
    throw std::runtime_error("playwright test bridge: no action found matching '" +
                             object_name_or_text + "'");
}

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

void send_mouse_press(SketcherWidget& sketcher, const std::string& object_name)
{
    // Squish send_event() delivers a Qt mouse press to the mapped object. Use
    // the same mechanism only when a QWidgetAction is not browser-rendered.
    auto& widget = require_visible_widget(sketcher, object_name);
    const QPointF point(widget.rect().center());
    QMouseEvent event(QEvent::MouseButtonPress, point, point,
                      widget.mapToGlobal(point.toPoint()), Qt::LeftButton,
                      Qt::LeftButton, Qt::NoModifier);
    QApplication::sendEvent(&widget, &event);
}

void close_active_popups()
{
    // Qt/WASM leaves QWidgetAction popup windows registered after the
    // Squish-style synthetic press. Close only the active popup stack, not
    // application windows, before returning to browser-visible input.
    while (auto* popup = QApplication::activePopupWidget()) {
        popup->close();
    }
}

std::string get_widget_rect(SketcherWidget& sketcher,
                            const std::string& object_name)
{
    auto* widget =
        find_visible_widget(sketcher, QString::fromStdString(object_name));
    if (!widget) {
        return "{}";
    }
    // Dialogs and menus are top-level Qt windows in WASM.  Map through global
    // coordinates so their browser-visible geometry is correct as well as
    // ordinary child-widget geometry.
    const QPoint top_left =
        sketcher.mapFromGlobal(widget->mapToGlobal(QPoint(0, 0)));
    return rect_to_json(top_left, widget->rect());
}

std::string get_widget_state(SketcherWidget& sketcher,
                             const std::string& object_name)
{
    auto& widget = require_widget(sketcher, object_name);
    QJsonObject result;
    result["objectName"] = widget.objectName();
    result["className"] = widget.metaObject()->className();
    result["visible"] = widget.isVisible();
    result["enabled"] = widget.isEnabled();
    result["styleSheet"] = widget.styleSheet();
    // Tool-selection suites assert the visible Qt tooltip text.  This is
    // read-only test observation and stays behind the Playwright bridge flag.
    result["toolTip"] = widget.toolTip();

    if (auto* button = qobject_cast<QAbstractButton*>(&widget)) {
        result["checked"] = button->isChecked();
        result["text"] = button->text();
    } else if (auto* line_edit = qobject_cast<QLineEdit*>(&widget)) {
        result["text"] = line_edit->text();
    } else if (auto* text_edit = qobject_cast<QPlainTextEdit*>(&widget)) {
        result["text"] = text_edit->toPlainText();
    } else if (auto* label = qobject_cast<QLabel*>(&widget)) {
        result["text"] = label->text();
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

void record_file_export(const std::string& filename,
                        const std::string& file_content)
{
#ifdef __EMSCRIPTEN__
    // QFileDialog::saveFileContent is implemented within Qt/WASM and does not
    // create a Chromium Download event.  Observe the already-produced bytes
    // here; the test itself still activates the visible Download button.
    QJsonObject result;
    result["filename"] = QString::fromStdString(filename);
    result["contentBase64"] = QString::fromLatin1(
        QByteArray::fromStdString(file_content).toBase64());
    const auto json = QJsonDocument(result).toJson(QJsonDocument::Compact);
    EM_ASM(
        {
            if (!window.__sketcherPlaywrightFileExports) {
                window.__sketcherPlaywrightFileExports = [];
            }
            window.__sketcherPlaywrightFileExports.push(
                JSON.parse(UTF8ToString($0)));
        },
        json.constData());
#else
    (void)filename;
    (void)file_content;
#endif
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

std::string get_rendered_atom_geometry(SketcherWidget& sketcher)
{
    QJsonArray result;
    for (auto* item : require_view(sketcher).scene()->items()) {
        auto* atom_item = qgraphicsitem_cast<AtomItem*>(item);
        if (!atom_item)
            continue;
        const auto* atom = atom_item->getAtom();
        const auto& point =
            atom->getOwningMol().getConformer().getAtomPos(atom->getIdx());
        QJsonObject entry;
        entry["index"] = static_cast<int>(atom->getIdx()) + 1;
        entry["x"] = point.x;
        entry["y"] = point.y;
        result.append(entry);
    }
    return QJsonDocument(result).toJson(QJsonDocument::Compact).toStdString();
}

std::string get_rendered_bond_geometry(SketcherWidget& sketcher)
{
    QJsonArray result;
    for (auto* item : require_view(sketcher).scene()->items()) {
        auto* bond_item = qgraphicsitem_cast<BondItem*>(item);
        if (!bond_item)
            continue;
        const auto* bond = bond_item->getBond();
        QJsonObject entry;
        entry["index"] = static_cast<int>(bond->getIdx()) + 1;
        entry["atom1"] = static_cast<int>(bond->getBeginAtomIdx()) + 1;
        entry["atom2"] = static_cast<int>(bond->getEndAtomIdx()) + 1;
        result.append(entry);
    }
    return QJsonDocument(result).toJson(QJsonDocument::Compact).toStdString();
}

std::string visible_widget_names(SketcherWidget& sketcher)
{
    QJsonArray result;
    for (auto* widget : sketcher.findChildren<QWidget*>()) {
        if (widget->isVisible() && !widget->objectName().isEmpty()) {
            result.append(widget->objectName());
        }
    }
    for (auto* top_level : QApplication::topLevelWidgets()) {
        if (!top_level || !top_level->isVisible())
            continue;
        if (!top_level->objectName().isEmpty())
            result.append(top_level->objectName());
        for (auto* widget : top_level->findChildren<QWidget*>()) {
            if (widget->isVisible() && !widget->objectName().isEmpty())
                result.append(widget->objectName());
        }
    }
    return QJsonDocument(result).toJson(QJsonDocument::Compact).toStdString();
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
