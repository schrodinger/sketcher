/* -------------------------------------------------------------------------
 * Playwright-only test bridge for the WebAssembly Sketcher application.
 *
 * Qt/WASM paints the entire application into a single canvas element, so
 * Playwright cannot see individual widgets or scene items in the DOM. The
 * bridge exists to close exactly that gap and no more: sketcher_get_rect()
 * resolves a selector to canvas coordinates, and the test then drives the
 * application with real Playwright mouse and keyboard events. Anything a test
 * can do with coordinates plus Playwright belongs in the JavaScript helpers,
 * not here.
 *
 * A Qt::Popup widget is painted into its own canvas element, but that element
 * is positioned over the page, so a rect mapped through global coordinates is
 * still a clickable page coordinate. Such a widget is reached with a real mouse
 * event like anything else.
 *
 * A context menu behaves the same way: the sketcher shows one with
 * QMenu::show() rather than exec(), so it does not run a nested event loop and
 * its rows are located with a "menu:" selector and clicked for real.
 *
 * A QToolButton's menu is the one exception, and sketcher_activate_action()
 * exists for it alone: Qt does run a nested event loop while such a menu is
 * open, which under Asyncify suspends the WebAssembly stack and stops any call
 * into the module from completing until the menu closes. See that function.
 *
 * Everything here is self-contained: the functions are registered with
 * JavaScript by the EMSCRIPTEN_BINDINGS block at the bottom, so no other
 * translation unit refers to them. Production UI code neither calls nor depends
 * on this interface. The bound names are underscore-prefixed to mark them
 * test-only.
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#ifdef __EMSCRIPTEN__
#include <emscripten/bind.h>
#endif

#include <stdexcept>
#include <string>
#include <utility>

#include <QAbstractButton>
#include <QAction>
#include <QApplication>
#include <QGraphicsItem>
#include <QGraphicsScene>
#include <QGraphicsView>
#include <QJsonDocument>
#include <QJsonObject>
#include <QMenu>
#include <QPoint>
#include <QRect>
#include <QSize>
#include <QString>
#include <QWidget>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>

#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"
#include "schrodinger/sketcher/molviewer/abstract_bond_or_connector_item.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/sketcher/widget/tool_button_with_popup.h"
#include "sketcher_instance.h"

using schrodinger::sketcher::AbstractAtomOrMonomerItem;
using schrodinger::sketcher::AbstractBondOrConnectorItem;
using schrodinger::sketcher::SketcherWidget;
using schrodinger::sketcher::ToolButtonWithPopup;

namespace
{

std::string to_json(const QJsonObject& object)
{
    return QJsonDocument(object).toJson(QJsonDocument::Compact).toStdString();
}

QJsonObject rect_json(const QPoint& top_left, const QSize& size,
                      const bool enabled)
{
    QJsonObject result;
    result["x"] = top_left.x();
    result["y"] = top_left.y();
    result["width"] = size.width();
    result["height"] = size.height();
    result["enabled"] = enabled;
    return result;
}

std::string rect_to_json(const QPoint& top_left, const QSize& size,
                         const bool enabled)
{
    return to_json(rect_json(top_left, size, enabled));
}

/**
 * Split a selector such as "widget:na_u_btn" or "atom:2" into its kind and
 * value. Throws if there is no ':' separator, which indicates a test bug.
 */
std::pair<std::string, std::string> split_selector(const std::string& selector)
{
    const auto separator = selector.find(':');
    if (separator == std::string::npos) {
        throw std::runtime_error(
            "playwright test bridge: malformed selector '" + selector +
            "' (expected \"<kind>:<value>\")");
    }
    return {selector.substr(0, separator), selector.substr(separator + 1)};
}

/**
 * QAction text may contain an ampersand marking its mnemonic (e.g. "Modify
 * &All"), but tests refer to actions by the name the user sees.
 */
QString without_mnemonic(QString text)
{
    text.remove('&');
    return text.simplified();
}

/**
 * Find a visible child widget by objectName. Several widgets may share a name,
 * so this returns the first visible one, or nullptr if none is visible.
 */
QWidget* find_visible_widget(SketcherWidget& sketcher, const QString& name)
{
    for (auto* widget : sketcher.findChildren<QWidget*>(name)) {
        if (widget->isVisible()) {
            return widget;
        }
    }
    return nullptr;
}

/**
 * Find the button whose popup contains the named widget, or "" if the widget
 * isn't in a popup at all.
 *
 * Many tools live in a popup that only appears while its button is pressed and
 * held, so a test that wants such a tool has to know which button to hold.
 * Deriving that from the widget tree keeps tests from carrying a hand-written
 * map of every tool to its owning button.
 */
std::string popup_owner(SketcherWidget& sketcher, const std::string& name)
{
    const auto wanted = QString::fromStdString(name);
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

/**
 * Map a point from some widget's coordinates into the sketcher's.
 *
 * This deliberately goes through global coordinates rather than using
 * QWidget::mapTo(). A Qt::Popup (or a dialog) is a separate top-level window,
 * and mapTo() only walks the widget hierarchy — for a window it accumulates
 * screen coordinates partway through and returns nonsense. Going via global
 * coordinates is correct for both cases and matches what Qt itself does for
 * independent top-level widgets.
 */
QPoint map_to_sketcher(const QWidget& widget, const SketcherWidget& sketcher,
                       const QPoint& point)
{
    return sketcher.mapFromGlobal(widget.mapToGlobal(point));
}

QGraphicsView& require_view(SketcherWidget& sketcher)
{
    auto* view = sketcher.findChild<QGraphicsView*>("view");
    if (view == nullptr) {
        throw std::runtime_error(
            "playwright test bridge: no QGraphicsView named 'view'");
    }
    return *view;
}

/**
 * Find the visible Scene item for the atom or bond at the given model index, or
 * nullptr if there isn't one.
 *
 * The casts are dynamic_casts to the abstract base classes rather than
 * qgraphicsitem_casts to AtomItem/BondItem, because qgraphicsitem_cast matches
 * type() exactly and so would skip monomers and monomer connectors: those are
 * siblings of AtomItem and BondItem, not subclasses.
 */
QGraphicsItem* find_visible_item(const QGraphicsView& view, const bool is_atom,
                                 const int index)
{
    for (auto* item : view.scene()->items()) {
        if (!item->isVisible()) {
            continue;
        }
        if (is_atom) {
            auto* atom_item = dynamic_cast<AbstractAtomOrMonomerItem*>(item);
            if (atom_item != nullptr &&
                static_cast<int>(atom_item->getAtom()->getIdx()) == index) {
                return atom_item;
            }
        } else {
            auto* bond_item = dynamic_cast<AbstractBondOrConnectorItem*>(item);
            if (bond_item != nullptr &&
                static_cast<int>(bond_item->getBond()->getIdx()) == index) {
                return bond_item;
            }
        }
    }
    return nullptr;
}

/**
 * Resolve a "widget:<objectName>" or "state:<objectName>" selector, or "{}" if
 * nothing matches.
 *
 * "widget" only matches a visible widget, since that is what a test can click.
 * "state" matches whether or not the widget is showing and reports "visible",
 * which is how a test asserts that a shortcut selected a tool that lives in a
 * closed popup. A button also reports "checked" and "text".
 */
std::string widget_rect(SketcherWidget& sketcher, const std::string& name,
                        const bool visible_only)
{
    const auto object_name = QString::fromStdString(name);
    auto* widget = find_visible_widget(sketcher, object_name);
    if (widget == nullptr) {
        if (visible_only) {
            return "{}";
        }
        widget = sketcher.findChild<QWidget*>(object_name);
        if (widget == nullptr) {
            return "{}";
        }
    }
    auto result = rect_json(map_to_sketcher(*widget, sketcher, QPoint(0, 0)),
                            widget->size(), widget->isEnabled());
    if (!visible_only) {
        result["visible"] = widget->isVisible();
    }
    if (auto* button = qobject_cast<QAbstractButton*>(widget)) {
        result["checked"] = button->isChecked();
        result["text"] = without_mnemonic(button->text());
        result["toolTip"] = button->toolTip();
    }
    return to_json(result);
}

/**
 * Resolve an "atom:<index>" or "bond:<index>" selector, or "{}" if no visible
 * item matches. Scene items are QGraphicsItems rather than QWidgets, so they
 * can't be found by objectName; their bounding rect is mapped through the View
 * transform instead.
 */
std::string item_rect(SketcherWidget& sketcher, const bool is_atom,
                      const std::string& value)
{
    bool is_number = false;
    const int index = QString::fromStdString(value).toInt(&is_number);
    if (!is_number) {
        throw std::runtime_error("playwright test bridge: '" + value +
                                 "' is not a valid item index");
    }
    auto& view = require_view(sketcher);
    auto* item = find_visible_item(view, is_atom, index);
    if (item == nullptr) {
        return "{}";
    }
    const QRect viewport_rect =
        view.mapFromScene(item->sceneBoundingRect()).boundingRect();
    return rect_to_json(
        map_to_sketcher(*view.viewport(), sketcher, viewport_rect.topLeft()),
        viewport_rect.size(), item->isEnabled());
}

/**
 * Resolve a row of one open QMenu, or "" if this menu has no matching row.
 *
 * Reading action geometry is only safe once the popup is actually visible; Qt
 * has not laid the rows out before that, and on WASM touching a menu's
 * internals mid-show aborts the runtime.
 */
std::string action_rect_in_menu(const QMenu& menu,
                                const SketcherWidget& sketcher,
                                const QString& name)
{
    if (!menu.isVisible()) {
        return {};
    }
    const QString wanted = without_mnemonic(name);
    for (auto* action : menu.actions()) {
        if (action->objectName() != name &&
            without_mnemonic(action->text()) != wanted) {
            continue;
        }
        const QRect rect = menu.actionGeometry(action);
        if (rect.isEmpty()) {
            continue;
        }
        return rect_to_json(map_to_sketcher(menu, sketcher, rect.topLeft()),
                            rect.size(), action->isEnabled());
    }
    return {};
}

/**
 * Resolve a "menu:<objectName or text>" selector against whichever menus are
 * open, or "{}" if none of them has a matching row.
 *
 * The active popup is tried first. With a submenu open Qt reports the deepest
 * one, which is both the menu a human pointer event would reach next and the
 * disambiguation this needs, since the same label can appear at more than one
 * level of a nested menu.
 */
std::string menu_rect(SketcherWidget& sketcher, const std::string& value)
{
    const auto name = QString::fromStdString(value);
    if (auto* active =
            qobject_cast<QMenu*>(QApplication::activePopupWidget())) {
        if (const auto result = action_rect_in_menu(*active, sketcher, name);
            !result.empty()) {
            return result;
        }
    }
    for (auto* menu : sketcher.findChildren<QMenu*>()) {
        if (const auto result = action_rect_in_menu(*menu, sketcher, name);
            !result.empty()) {
            return result;
        }
    }
    return "{}";
}

/**
 * Trigger the QAction matching name_or_text, or return false if there is none.
 *
 * Actions are matched on objectName first and then on display text, since most
 * menu actions are created without an objectName.
 */
bool try_activate_action(SketcherWidget& sketcher, const QString& name)
{
    const QString wanted = without_mnemonic(name);
    for (auto* action : sketcher.findChildren<QAction*>()) {
        if (action->objectName() != name &&
            without_mnemonic(action->text()) != wanted) {
            continue;
        }
        if (!action->isEnabled()) {
            throw std::runtime_error("playwright test bridge: action '" +
                                     name.toStdString() + "' is disabled");
        }
        action->trigger();
        return true;
    }
    return false;
}

} // namespace

// The functions below deliberately have external linkage even though nothing
// else refers to them: the EMSCRIPTEN_BINDINGS block is compiled out on desktop
// builds, and file-local functions would then trip -Wunused-function under
// -Werror. Compiling them everywhere means the non-WASM CI jobs still catch
// errors here.

/**
 * Resolve an object selector to its position and size on the canvas, as
 * {"x":…, "y":…, "width":…, "height":…, "enabled":…}. Coordinates are relative
 * to the sketcher widget's top-left corner, which is also the top-left of the
 * WASM canvas, so tests can use them directly as page coordinates.
 *
 * Supported selectors:
 *
 *   "widget:<objectName>"  a visible QWidget, by its Qt objectName
 *   "state:<objectName>"   the same, but matching a hidden widget too
 *   "atom:<index>"         an atom of the current molecule, by model index
 *   "bond:<index>"         a bond of the current molecule, by model index
 *   "menu:<name or text>"  a row of an open context menu, by objectName or text
 *
 * Use "widget:" to find something to click, since only a visible widget can be
 * clicked. Use "state:" to ask whether a control is checked or enabled when it
 * may be out of sight, such as a tool inside a closed popup; it adds "visible"
 * to the result. A button of either kind also reports "checked", "text", and
 * "toolTip".
 *
 * A "menu:" selector only resolves while the menu is on screen, so open the
 * menu first. It does not reach a QToolButton's menu, which cannot be open and
 * queried at the same time; use sketcher_activate_action() for those rows.
 *
 * Monomers are addressed as "atom" and monomer connectors as "bond", since the
 * model stores them as RDKit atoms and bonds. Every atom has a non-empty
 * bounding rect even when its label isn't painted (an unlabeled carbon, say),
 * because the rect is built from the predictive highlighting path, so the
 * returned rect is always a usable click target.
 *
 * Model indices are indices into the molecule returned by
 * SketcherWidget::getRDKitMolecule(). That molecule is a plain copy of the one
 * the scene items are built from, so the indices agree; note that this is the
 * order the structure was built in, not the canonical order a format like
 * SMILES may renumber to on export.
 *
 * Returns "{}" when nothing visible matches, since there is no coordinate a
 * test could click. "enabled" reports whether Qt would accept a click there;
 * Playwright applies that check automatically to real DOM elements but has no
 * way to do so for something painted into a canvas.
 *
 * Throws std::runtime_error for a malformed selector or an unrecognized kind,
 * which indicates the test needs to be updated.
 */
std::string sketcher_get_rect(const std::string& selector)
{
    const auto [kind, value] = split_selector(selector);
    auto& sketcher = get_sketcher_instance();
    if (kind == "widget" || kind == "state") {
        return widget_rect(sketcher, value, kind == "widget");
    }
    if (kind == "atom" || kind == "bond") {
        return item_rect(sketcher, kind == "atom", value);
    }
    if (kind == "menu") {
        return menu_rect(sketcher, value);
    }
    throw std::runtime_error(
        "playwright test bridge: unrecognized selector kind '" + kind +
        "' (expected 'widget', 'state', 'atom', 'bond', or 'menu')");
}

/**
 * Trigger a menu action by objectName or visible text, without opening the menu
 * it belongs to. Mnemonic ampersands are ignored when comparing text.
 *
 * This exists because a menu row cannot be clicked the way every other control
 * can. The menu belongs to a QToolButton, and Qt runs a nested event loop for
 * as long as that menu is open. Qt/WASM is built with Asyncify, so that nested
 * loop leaves the WebAssembly stack suspended, and no further call into the
 * module completes until the menu closes — sketcher_get_rect() included. A test
 * therefore cannot ask where a menu row is while the menu is showing, which
 * leaves triggering the action directly as the only way to reach it.
 *
 * Note what this gives up: the menu never opens, so nothing here covers the
 * menu's own layout or hit-testing, only the action's effect. Controls outside
 * a menu — including those in Qt::Popup widgets, which do not run a nested
 * event loop — should still be clicked with real mouse events.
 *
 * Throws std::runtime_error if nothing matches, or if the match is disabled.
 * Triggering skips the enabled check a real mouse event goes through, so
 * refusing here keeps a test from passing against a row the user could not have
 * activated.
 */
void sketcher_activate_action(const std::string& name_or_text)
{
    auto& sketcher = get_sketcher_instance();
    if (!try_activate_action(sketcher, QString::fromStdString(name_or_text))) {
        throw std::runtime_error(
            "playwright test bridge: no action found matching '" +
            name_or_text + "'");
    }
}

/**
 * Return the objectName of the button whose popup holds the named widget, or an
 * empty string if it isn't in a popup.
 *
 * A test uses this to find the button it has to press and hold before a tool
 * becomes visible and clickable.
 */
std::string sketcher_get_popup_owner(const std::string& name)
{
    return popup_owner(get_sketcher_instance(), name);
}

#ifdef __EMSCRIPTEN__
// A second bindings block alongside the one in main.cpp; embind allows any
// number of them as long as each has a distinct name. This object file is
// linked directly into the executable rather than through a static library, so
// the registrations always run.
EMSCRIPTEN_BINDINGS(sketcher_playwright_test_bridge)
{
    emscripten::function("_sketcher_get_rect", &sketcher_get_rect);
    emscripten::function("_sketcher_activate_action",
                         &sketcher_activate_action);
    emscripten::function("_sketcher_get_popup_owner",
                         &sketcher_get_popup_owner);
}
#endif
