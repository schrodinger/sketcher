#include "schrodinger/sketcher/widget/multi_select_combo_box.h"

#include <algorithm>

#include <QAbstractItemView>
#include <QApplication>
#include <QFontMetrics>
#include <QHoverEvent>
#include <QKeyEvent>
#include <QListView>
#include <QMouseEvent>
#include <QPalette>
#include <QSignalBlocker>
#include <QStandardItemModel>
#include <QStringList>
#include <QStyle>
#include <QStyleOptionComboBox>
#include <QStylePainter>

namespace schrodinger
{
namespace sketcher
{

static const QString ALL_CHECKED_TEXT = "Any";
static const QString NONE_CHECKED_TEXT = "None";
static const QString DISPLAY_TEXT_SEPARATOR = ", ";
// The empty space around each item in the dropdown, in pixels
static const int ITEM_SPACING = 2;

/**
 * @return Whether the key toggles the current item of the dropdown
 */
static bool is_toggle_key(int key)
{
    return key == Qt::Key_Space || key == Qt::Key_Enter ||
           key == Qt::Key_Return;
}

/**
 * @return Whether the key moves the dropdown's current item
 */
static bool is_navigation_key(int key)
{
    return key == Qt::Key_Up || key == Qt::Key_Down || key == Qt::Key_Home ||
           key == Qt::Key_End || key == Qt::Key_PageUp ||
           key == Qt::Key_PageDown;
}

/**
 * @return Whether the key event has no modifiers, ignoring the keypad modifier
 * so that the keypad's Enter and arrow keys behave like the main ones
 */
static bool has_no_modifiers(const QKeyEvent* event)
{
    return (event->modifiers() & ~Qt::KeypadModifier) == Qt::NoModifier;
}

MultiSelectComboBox::MultiSelectComboBox(QWidget* parent) : QComboBox(parent)
{
    m_model = new QStandardItemModel(this);
    setModel(m_model);
    connect(m_model, &QStandardItemModel::rowsInserted, this,
            &MultiSelectComboBox::onRowsInserted);
    connect(m_model, &QStandardItemModel::dataChanged, this,
            &MultiSelectComboBox::onDataChanged);
    connect(m_model, &QStandardItemModel::rowsRemoved, this,
            &MultiSelectComboBox::checkedItemsChanged);
    // The display text depends on the check states, which don't change the
    // current item, so repaint whenever they change
    connect(this, &MultiSelectComboBox::checkedItemsChanged, this,
            qOverload<>(&MultiSelectComboBox::update));
    // Calling view() creates the dropdown, which installs its own event
    // filters on the view and its viewport. Event filters installed later run
    // first, so these filters see each event before the dropdown can close
    // itself in response to it.
    view()->installEventFilter(this);
    view()->viewport()->installEventFilter(this);
    // The item checkboxes fill most of each item's height, so space out the
    // items to keep the checkboxes away from each other and from the edges of
    // the dropdown. QComboBox accounts for the spacing when it sizes the
    // dropdown.
    if (auto* list_view = qobject_cast<QListView*>(view())) {
        list_view->setSpacing(ITEM_SPACING);
    }
}

void MultiSelectComboBox::onRowsInserted(const QModelIndex& parent, int first,
                                         int last)
{
    if (parent.isValid()) {
        return;
    }
    QSignalBlocker blocker(m_model);
    for (int row = first; row <= last; ++row) {
        auto* item = m_model->item(row);
        // Items are selectable so that the view's selection can highlight the
        // current item, which is the item that keyboard input applies to.
        // The selection doesn't affect which items are checked.
        item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsSelectable |
                       Qt::ItemIsUserCheckable);
        item->setCheckState(Qt::Unchecked);
    }
}

void MultiSelectComboBox::onDataChanged(const QModelIndex&, const QModelIndex&,
                                        const QList<int>& roles)
{
    if (!roles.isEmpty() && !roles.contains(Qt::CheckStateRole)) {
        return;
    }
    if (m_suppress_changed_signal) {
        m_changed_while_suppressed = true;
    } else {
        emit checkedItemsChanged();
    }
}

void MultiSelectComboBox::setCheckedData(const QList<QVariant>& data)
{
    m_suppress_changed_signal = true;
    m_changed_while_suppressed = false;
    for (int row = 0; row < count(); ++row) {
        setItemChecked(row, data.contains(itemData(row)));
    }
    m_suppress_changed_signal = false;
    if (m_changed_while_suppressed) {
        emit checkedItemsChanged();
    }
}

QList<QVariant> MultiSelectComboBox::checkedData() const
{
    QList<QVariant> data;
    for (int row = 0; row < count(); ++row) {
        if (isItemChecked(row)) {
            data.append(itemData(row));
        }
    }
    return data;
}

void MultiSelectComboBox::setItemChecked(int index, bool checked)
{
    auto* item = m_model->item(index);
    if (item == nullptr) {
        return;
    }
    const auto state = checked ? Qt::Checked : Qt::Unchecked;
    if (item->checkState() != state) {
        item->setCheckState(state);
    }
}

bool MultiSelectComboBox::isItemChecked(int index) const
{
    auto* item = m_model->item(index);
    return item != nullptr && item->checkState() == Qt::Checked;
}

QString MultiSelectComboBox::getDisplayText() const
{
    QStringList checked_texts;
    for (int row = 0; row < count(); ++row) {
        if (isItemChecked(row)) {
            checked_texts.append(itemText(row));
        }
    }
    if (checked_texts.isEmpty()) {
        return NONE_CHECKED_TEXT;
    } else if (checked_texts.size() == count()) {
        return ALL_CHECKED_TEXT;
    }
    return checked_texts.join(DISPLAY_TEXT_SEPARATOR);
}

bool MultiSelectComboBox::toggleItemByUser(int row)
{
    auto* item = m_model->item(row);
    if (item == nullptr || !item->isEnabled()) {
        return false;
    }
    setItemChecked(row, !isItemChecked(row));
    return true;
}

bool MultiSelectComboBox::isCurrentItemHighlighted() const
{
    return view()->selectionModel()->isSelected(view()->currentIndex());
}

void MultiSelectComboBox::highlightCurrentItem()
{
    view()->selectionModel()->select(view()->currentIndex(),
                                     QItemSelectionModel::ClearAndSelect);
}

bool MultiSelectComboBox::handleDropdownKeyPress(QKeyEvent* event)
{
    const int key = event->key();
    if (is_navigation_key(key)) {
        // The mouse may still be hovering over an item that some styles
        // highlight, so clear that highlight to leave only the keyboard's
        // highlight
        clearHoverHighlight();
        if ((key == Qt::Key_Up || key == Qt::Key_Down) &&
            view()->currentIndex().isValid() && !isCurrentItemHighlighted()) {
            // The user couldn't see which item was current, so reveal it
            // rather than moving away from it
            highlightCurrentItem();
            return true;
        }
        // Let the view move the current item, which also selects it
        return false;
    }
    if (!is_toggle_key(key)) {
        return false;
    }
    if (event->isAutoRepeat()) {
        return true;
    }
    const bool toggled = toggleItemByUser(view()->currentIndex().row());
    if (key == Qt::Key_Space) {
        highlightCurrentItem();
    } else if (toggled) {
        // Enter and Return also close the dropdown, so the user can quickly
        // toggle a single item
        hidePopup();
    }
    return true;
}

bool MultiSelectComboBox::eventFilter(QObject* watched, QEvent* event)
{
    auto* item_view = view();
    if (watched == item_view->viewport()) {
        switch (event->type()) {
            case QEvent::MouseButtonPress:
            case QEvent::MouseButtonDblClick: {
                auto* mouse_event = static_cast<QMouseEvent*>(event);
                m_pressed_row =
                    item_view->indexAt(mouse_event->position().toPoint()).row();
                // Let a single press through so that the view can make the
                // pressed row current, but swallow a double-click so that it
                // can't toggle the item on its own
                return event->type() == QEvent::MouseButtonDblClick;
            }
            case QEvent::MouseButtonRelease: {
                auto* mouse_event = static_cast<QMouseEvent*>(event);
                const int row =
                    item_view->indexAt(mouse_event->position().toPoint()).row();
                if (row >= 0 && row == m_pressed_row) {
                    toggleItemByUser(row);
                }
                m_pressed_row = -1;
                // Swallow the release so that the dropdown stays open and the
                // view's delegate doesn't also toggle the item
                return true;
            }
            default:
                break;
        }
    } else if (watched == item_view) {
        if (event->type() == QEvent::ShortcutOverride ||
            event->type() == QEvent::KeyPress) {
            auto* key_event = static_cast<QKeyEvent*>(event);
            if (has_no_modifiers(key_event)) {
                if (event->type() == QEvent::KeyPress) {
                    return handleDropdownKeyPress(key_event);
                } else if (is_toggle_key(key_event->key())) {
                    // Ask for the key press rather than letting the dropdown
                    // close itself in response to Enter or Return without
                    // toggling the item
                    event->accept();
                    return true;
                }
            }
        }
    }
    return QComboBox::eventFilter(watched, event);
}

void MultiSelectComboBox::clearHoverHighlight()
{
    QHoverEvent leave(QEvent::HoverLeave, QPointF(-1, -1), QPointF(-1, -1),
                      QPointF(-1, -1));
    QApplication::sendEvent(view()->viewport(), &leave);
}

void MultiSelectComboBox::updateDropdownBackground()
{
    QStyleOptionComboBox option;
    initStyleOption(&option);
    auto* viewport = view()->viewport();
    if (!style()->styleHint(QStyle::SH_ComboBox_Popup, &option, this)) {
        // The items have no background of their own, so the view's base color
        // already shows behind them
        viewport->setBackgroundRole(QPalette::Base);
        return;
    }
    // Use the same palette that Qt's menu-style item delegate
    // (QComboMenuDelegate) uses to fill each item
    const auto item_palette =
        view()->palette().resolve(QApplication::palette("QMenu"));
    auto viewport_palette = viewport->palette();
    for (auto group :
         {QPalette::Active, QPalette::Inactive, QPalette::Disabled}) {
        viewport_palette.setBrush(group, QPalette::Window,
                                  item_palette.brush(group, QPalette::Window));
    }
    viewport->setPalette(viewport_palette);
    viewport->setBackgroundRole(QPalette::Window);
}

void MultiSelectComboBox::showPopup()
{
    // The style or palette may have changed since the dropdown was last shown
    updateDropdownBackground();
    QComboBox::showPopup();
    clearHoverHighlight();
    // Keep the current item, which keyboard input applies to, but don't
    // highlight it until the user hovers over an item or presses a key
    view()->selectionModel()->clearSelection();
}

void MultiSelectComboBox::hidePopup()
{
    QComboBox::hidePopup();
    clearHoverHighlight();
    m_pressed_row = -1;
}

QSize MultiSelectComboBox::sizeHint() const
{
    auto hint = QComboBox::sizeHint();
    // The widest text that can be displayed has every item checked but one
    const auto metrics = fontMetrics();
    int text_width = std::max(metrics.horizontalAdvance(ALL_CHECKED_TEXT),
                              metrics.horizontalAdvance(NONE_CHECKED_TEXT));
    for (int skipped = 0; skipped < count(); ++skipped) {
        QStringList texts;
        for (int row = 0; row < count(); ++row) {
            if (row != skipped) {
                texts.append(itemText(row));
            }
        }
        text_width = std::max(text_width, metrics.horizontalAdvance(texts.join(
                                              DISPLAY_TEXT_SEPARATOR)));
    }
    QStyleOptionComboBox option;
    initStyleOption(&option);
    const auto size_from_contents = style()->sizeFromContents(
        QStyle::CT_ComboBox, &option, QSize(text_width, hint.height()), this);
    hint.setWidth(std::max(hint.width(), size_from_contents.width()));
    return hint;
}

void MultiSelectComboBox::paintEvent(QPaintEvent*)
{
    QStylePainter painter(this);
    painter.setPen(palette().color(QPalette::Text));
    QStyleOptionComboBox option;
    initStyleOption(&option);
    option.currentIcon = QIcon();
    const auto text_rect = style()->subControlRect(
        QStyle::CC_ComboBox, &option, QStyle::SC_ComboBoxEditField, this);
    option.currentText = fontMetrics().elidedText(
        getDisplayText(), Qt::ElideRight, text_rect.width());
    painter.drawComplexControl(QStyle::CC_ComboBox, option);
    painter.drawControl(QStyle::CE_ComboBoxLabel, option);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/multi_select_combo_box.moc"
