#include "schrodinger/sketcher/menu/sketcher_top_bar_menus.h"

#include <vector>

#include <QKeySequence>
#include <QString>

#include "schrodinger/sketcher/dialog/about_2d_sketcher.h"
#include "schrodinger/sketcher/dialog/sketcher_welcome_dialog.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

namespace
{

/**
 * Convenience function for adding an action to a menu with a visible shortcut.
 *
 * @param menu The menu to which the action should be added
 * @param text The text to display on the new action
 * @param shortcut The shortcut sequence
 * @return A pointer to the new action
 */
void add_action_with_shortcut(QMenu* menu, QAction* action,
                              const QKeySequence& shortcut)
{
    menu->addAction(action);
    action->setShortcut(shortcut);
    action->setShortcutVisibleInContextMenu(true);
}

QAction* add_action_with_shortcut(QMenu* menu, const QString& text,
                                  const QKeySequence& shortcut)
{
    auto action = new QAction(text, menu);
    add_action_with_shortcut(menu, action, shortcut);
    return action;
}
} // namespace

namespace schrodinger
{
namespace sketcher
{

ImportMenu::ImportMenu(QWidget* parent) : QMenu(parent)
{
    m_import_from_file_act = addAction("Import from File...");
    m_paste_in_text_act = addAction("Paste in Text...");
    m_replace_content_act = addAction("Replace Current Content");
    m_replace_content_act->setCheckable(true);
}

ExportMenu::ExportMenu(QWidget* parent) : QMenu(parent)
{
    m_save_image_act = addAction("Save Image...");
    m_export_to_file_act = addAction("Export to File...");
}

MoreActionsMenu::MoreActionsMenu(SketcherModel* model, QWidget* parent) :
    QMenu(parent)
{
    m_modify_structure_menu = addMenu("Modify All");
    addSeparator();
    m_undo_act = add_action_with_shortcut(this, "Undo", tr("Ctrl+Z"));
    m_redo_act = add_action_with_shortcut(this, "Redo", tr("Shift+Ctrl+Z"));
    addSeparator();
    m_select_all_act =
        add_action_with_shortcut(this, "Select All", tr("Ctrl+A"));
    m_clear_selection_act =
        add_action_with_shortcut(this, "Clear Selection", tr("Ctrl+D"));
    m_expand_selection_menu = addMenu("Expand Selection");
    m_invert_selection_act =
        add_action_with_shortcut(this, "Invert Selection", tr("Ctrl+I"));
    addSeparator();

    m_cut_copy_manager = new CutCopyActionManager(this);
    add_action_with_shortcut(this, m_cut_copy_manager->m_cut_action,
                             tr("Ctrl+X"));
    add_action_with_shortcut(this, m_cut_copy_manager->m_copy_action,
                             tr("Ctrl+C"));
    addMenu(m_cut_copy_manager->m_copy_as_menu);
    m_paste_act = add_action_with_shortcut(this, "Paste", tr("Ctrl+V"));
    addSeparator();

    m_add_custom_fragment_act = addAction("Add Custom Fragment...");
    m_fit_to_screen_act =
        add_action_with_shortcut(this, "Fit to Screen", tr("Ctrl+F"));

    // "Modify All" submenu
    m_flip_horizontal_act =
        m_modify_structure_menu->addAction("Flip Horizontal");
    m_flip_vertical_act = m_modify_structure_menu->addAction("Flip Vertical");
    m_modify_structure_menu->addSeparator();
    m_aromatize_act = m_modify_structure_menu->addAction("Aromatize");
    m_kekulize_act = m_modify_structure_menu->addAction("Kekulize");
    m_modify_structure_menu->addSeparator();
    m_add_explicit_hydrogens_act =
        m_modify_structure_menu->addAction("Add Explicit Hydrogens");
    m_remove_explicit_hydrogens_act =
        m_modify_structure_menu->addAction("Remove Explicit Hydrogens");
}

void MoreActionsMenu::setModel(SketcherModel* model)
{
    m_cut_copy_manager->setModel(model);
}

ConfigureViewMenu::ConfigureViewMenu(QWidget* parent) : QMenu(parent)
{
    m_valence_errors_act = addAction("Valence Errors");
    m_heteroatom_colors_act = addAction("Heteroatom Colors");
    m_stereo_labels_act = addAction("Stereo Labels");
    m_implicit_hydrogens_act = addAction("Implicit Hydrogens");
    addSeparator();
    m_preferences_act = addAction("Preferences...");

    // Assign checkable actions
    std::vector<QAction*> checkable_actions = {
        m_valence_errors_act,
        m_heteroatom_colors_act,
        m_stereo_labels_act,
        m_implicit_hydrogens_act,
    };
    for (auto& action : checkable_actions) {
        action->setCheckable(true);
    }
}

HelpMenu::HelpMenu(QWidget* parent) : QMenu(parent)
{
    m_help_act = addAction("Help...");
    addAction("Getting Started...", this, &HelpMenu::showWelcomeDialog);
    addAction("About 2D Sketcher...", this,
              &HelpMenu::showAbout2DSketcherDialog);
}

void HelpMenu::showWelcomeDialog()
{
    auto dialog = new SketcherWelcomeDialog(parentWidget());
    dialog->setDoNotShowCheckboxVisible(false);
    dialog->show();
}

void HelpMenu::showAbout2DSketcherDialog()
{
    auto dialog = new About2DSketcher(parentWidget());
    dialog->show();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/sketcher_top_bar_menus.moc"
