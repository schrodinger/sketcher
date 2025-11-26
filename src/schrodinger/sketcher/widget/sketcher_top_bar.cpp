#include "schrodinger/sketcher/widget/sketcher_top_bar.h"

#include <QButtonGroup>
#include <QDesktopServices>
#include <QFileDialog>
#include <QRegularExpression>
#include <QToolButton>
#include <QWidget>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/rdkit_extensions/file_stream.h"
#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/dialog/paste_in_text_dialog.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/menu/sketcher_top_bar_menus.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/ui/ui_sketcher_top_bar.h"
#include "schrodinger/sketcher/version.h"
#include "schrodinger/sketcher/sketcher_css_style.h"

namespace schrodinger
{
namespace sketcher
{

static const char* APPLY_TO_SELECTION_PROPERTY = "apply_to_selection";

SketcherTopBar::SketcherTopBar(QWidget* parent) : SketcherView(parent)
{
    ui.reset(new Ui::SketcherTopBar());
    ui->setupUi(this);

    initButtons();
    initMenus();
}

SketcherTopBar::~SketcherTopBar() = default;

void SketcherTopBar::initButtons()
{
    connect(ui->undo_btn, &QToolButton::clicked, this,
            &SketcherTopBar::undoRequested);
    connect(ui->redo_btn, &QToolButton::clicked, this,
            &SketcherTopBar::redoRequested);
    connect(ui->fit_btn, &QToolButton::clicked, this, [this]() {
        bool apply_to_selection =
            ui->fit_btn->property(APPLY_TO_SELECTION_PROPERTY).toBool();
        emit fitToScreenRequested(apply_to_selection);
    });
    connect(ui->cleanup_btn, &QToolButton::clicked, this, [this]() {
        bool apply_to_selection =
            ui->cleanup_btn->property(APPLY_TO_SELECTION_PROPERTY).toBool();
        emit cleanupRequested(apply_to_selection);
    });
    connect(ui->clear_btn, &QToolButton::clicked, this,
            &SketcherTopBar::clearSketcherRequested);
}

void SketcherTopBar::initMenus()
{
    // Set up "New Structures" menu
    m_import_menu = new ImportMenu(this);
    ui->import_btn->setMenu(m_import_menu);
    connect(m_import_menu->m_import_from_file_act, &QAction::triggered, this,
            &SketcherTopBar::onImportFromFileClicked);
    connect(m_import_menu->m_paste_in_text_act, &QAction::triggered, this,
            &SketcherTopBar::onPasteInTextClicked);

    // Set up "Export" menu
    m_export_menu = new ExportMenu(this);
    ui->export_btn->setMenu(m_export_menu);
    connect(m_export_menu->m_save_image_act, &QAction::triggered, this,
            &SketcherTopBar::saveImageRequested);
    connect(m_export_menu->m_export_to_file_act, &QAction::triggered, this,
            &SketcherTopBar::exportToFileRequested);

    // Set up "More Actions" menu
    m_more_actions_menu = new MoreActionsMenu(getModel(), this);
    ui->more_actions_btn->setMenu(m_more_actions_menu);
    connect(m_more_actions_menu->m_undo_act, &QAction::triggered, this,
            &SketcherTopBar::undoRequested);
    connect(m_more_actions_menu->m_redo_act, &QAction::triggered, this,
            &SketcherTopBar::redoRequested);
    connect(m_more_actions_menu->m_select_all_act, &QAction::triggered, this,
            &SketcherTopBar::selectAllRequested);
    connect(m_more_actions_menu->m_clear_selection_act, &QAction::triggered,
            this, &SketcherTopBar::clearSelectionRequested);
    connect(m_more_actions_menu->m_cut_copy_manager,
            &CutCopyActionManager::cutRequested, this,
            &SketcherTopBar::cutRequested);
    connect(m_more_actions_menu->m_cut_copy_manager,
            &CutCopyActionManager::copyRequested, this,
            &SketcherTopBar::copyRequested);
    connect(m_more_actions_menu->m_paste_act, &QAction::triggered, this,
            &SketcherTopBar::pasteRequested);
    connect(m_more_actions_menu->m_flip_horizontal_act, &QAction::triggered,
            this, &SketcherTopBar::flipHorizontalRequested);
    connect(m_more_actions_menu->m_flip_vertical_act, &QAction::triggered, this,
            &SketcherTopBar::flipVerticalRequested);
    connect(m_more_actions_menu->m_add_explicit_hydrogens_act,
            &QAction::triggered, this,
            &SketcherTopBar::addExplicitHydrogensRequested);
    connect(m_more_actions_menu->m_aromatize_act, &QAction::triggered, this,
            &SketcherTopBar::aromatizeRequested);
    connect(m_more_actions_menu->m_kekulize_act, &QAction::triggered, this,
            &SketcherTopBar::kekulizeRequested);
    connect(m_more_actions_menu->m_remove_explicit_hydrogens_act,
            &QAction::triggered, this,
            &SketcherTopBar::removeExplicitHydrogensRequested);
    connect(m_more_actions_menu->m_add_custom_fragment_act, &QAction::triggered,
            this, &SketcherTopBar::onAddCustomFragmentClicked);
    connect(m_more_actions_menu->m_fit_to_screen_act, &QAction::triggered, this,
            &SketcherTopBar::fitToScreenRequested);

    // Set up "Configure View" menu
    m_configure_view_menu = new ConfigureViewMenu(this);
    ui->configure_view_btn->setMenu(m_configure_view_menu);
    connect(m_configure_view_menu->m_preferences_act, &QAction::triggered, this,
            &SketcherTopBar::adjustRenderingSettingsRequested);

    // Set up "Help" menu
    m_help_menu = new HelpMenu(this);
    ui->help_btn->setMenu(m_help_menu);
    connect(m_help_menu->m_help_act, &QAction::triggered, this,
            &SketcherTopBar::onHelpClicked);
    connect(m_more_actions_menu->m_invert_selection_act, &QAction::triggered,
            this, &SketcherTopBar::invertSelectionRequested);

    // FIXME: For now, hide menu items that have not been implemented
    for (auto action : {
             // TODO: SKETCH-1013
             m_more_actions_menu->m_expand_selection_menu->menuAction(),
             // TODO: SKETCH-1015
             m_more_actions_menu->m_add_custom_fragment_act,
             // TODO: SKETCH-1017
             m_configure_view_menu->m_implicit_hydrogens_act,
         }) {
        action->setVisible(false);
    }
}

void SketcherTopBar::setModel(SketcherModel* model)
{
    SketcherView::setModel(model);
    m_more_actions_menu->setModel(model);
    connect(model, &SketcherModel::interactiveItemsChanged, this,
            &SketcherTopBar::updateWidgetsEnabled);
    connect(model, &SketcherModel::undoStackDataChanged, this,
            &SketcherTopBar::updateWidgetsEnabled);
}

void SketcherTopBar::updateWidgetsEnabled()
{
    auto model = getModel();
    if (model == nullptr) {
        return;
    }

    // Selection based items should only be enabled with an active selection
    auto has_selection = model->hasActiveSelection();
    for (const auto& act : {
             m_more_actions_menu->m_expand_selection_menu->menuAction(),
             m_more_actions_menu->m_invert_selection_act,
             m_more_actions_menu->m_clear_selection_act,
         }) {
        act->setEnabled(has_selection);
    }
    ui->fit_btn->setProperty(APPLY_TO_SELECTION_PROPERTY, has_selection);
    ui->fit_btn->setStyleSheet(has_selection ? SELECTION_ACTIVE_STYLE : "");
    ui->cleanup_btn->setProperty(APPLY_TO_SELECTION_PROPERTY, has_selection);
    ui->cleanup_btn->setStyleSheet(has_selection ? SELECTION_ACTIVE_STYLE : "");

    auto [can_undo, can_redo] = model->getUndoStackData();
    ui->undo_btn->setEnabled(can_undo);
    m_more_actions_menu->m_undo_act->setEnabled(can_undo);
    ui->redo_btn->setEnabled(can_redo);
    m_more_actions_menu->m_redo_act->setEnabled(can_redo);

    // Modes, buttons, and menu items that are not applicable with empty scene
    auto scene_has_contents = !model->sceneIsEmpty();
    for (auto btn :
         {ui->fit_btn, ui->cleanup_btn, ui->clear_btn, ui->export_btn}) {
        btn->setEnabled(scene_has_contents);
    }
    for (auto act : {m_more_actions_menu->m_modify_structure_menu->menuAction(),
                     m_more_actions_menu->m_fit_to_screen_act}) {
        act->setEnabled(scene_has_contents);
    }
    m_more_actions_menu->m_select_all_act->setEnabled(
        scene_has_contents && !model->allItemsSelected());

    ui->clear_btn->setDisabled(!scene_has_contents);
}

void SketcherTopBar::generatePackets()
{
    m_signal_packets.emplace_back(ModelKey::NEW_STRUCTURES_REPLACE_CONTENT,
                                  m_import_menu->m_replace_content_act);
    m_signal_packets.emplace_back(ModelKey::SHOW_VALENCE_ERRORS,
                                  m_configure_view_menu->m_valence_errors_act);
    m_signal_packets.emplace_back(
        ModelKey::COLOR_HETEROATOMS,
        m_configure_view_menu->m_heteroatom_colors_act);
    m_signal_packets.emplace_back(ModelKey::SHOW_STEREO_LABELS,
                                  m_configure_view_menu->m_stereo_labels_act);
    m_signal_packets.emplace_back(
        ModelKey::USE_IMPLICIT_HYDROGENS,
        m_configure_view_menu->m_implicit_hydrogens_act);

    for (auto& signal_packet : m_signal_packets) {
        auto setter = std::bind(&QAction::setChecked, signal_packet.action,
                                std::placeholders::_1);
        m_setter_packets.emplace_back(signal_packet.key, setter);
    }
}

void SketcherTopBar::onImportFromFileClicked()
{
    auto file_open_completed = [this](const auto& file_path,
                                      const auto& content) {
        if (file_path.isEmpty()) {
            /*
             * If the user cancels the file dialog, the file path will be empty.
             * In this case, we do not want to do anything. This is the
             * behavior recommended in the QFileDialog::getOpenFileContent
             * documentation. SKETCH-2239
             */
            return;
        }
        try {
            using namespace rdkit_extensions;

            auto path_string = file_path.toStdString();
            auto format = get_file_format(path_string);
            auto text = get_decompressed_string(content.toStdString());
            emit importTextRequested(text, format);
        } catch (const std::exception& exc) {
            show_error_dialog("File Error", exc.what(), this);
        }
    };

    QStringList filters;
    for (const auto& [_, label, extensions] : get_import_formats()) {
        filters.append(get_filter_name(label, extensions));
    }
    auto name_filter = filters.join(";;");

#if QT_6_5_5
    QFileDialog::getOpenFileContent(name_filter, file_open_completed, this);
#else
    QFileDialog::getOpenFileContent(name_filter, file_open_completed);
#endif
}

void SketcherTopBar::onPasteInTextClicked()
{
    auto model = getModel();
    if (model == nullptr) {
        throw std::runtime_error(
            "The dialog cannot be used if there is no model.");
    }

    auto dialog = new PasteInTextDialog(model, this);
    connect(dialog, &PasteInTextDialog::textAccepted, this,
            &SketcherTopBar::importTextRequested);
    dialog->show();
}

void SketcherTopBar::onAddCustomFragmentClicked()
{
}

bool SketcherTopBar::handleShortcutAction(const QKeySequence& key_seq)
{
    for (auto action : m_more_actions_menu->actions()) {
        if (action->shortcut() == key_seq) {
            action->trigger();
            return true;
        }
    }
    // Undo/Redo button do not have action, but have their shortcuts.
    for (auto button : {ui->redo_btn, ui->undo_btn}) {
        if (button->shortcut() == key_seq) {
            button->click();
            return true;
        }
    }
    return false;
}

void SketcherTopBar::onHelpClicked()
{
    auto documentation_url =
        QString("https://learn.schrodinger.com/public/2D-Sketcher/%1/Content/"
                "2d-sketcher/2d_sketcher_home.htm")
            .arg(std::string(SKETCHER_RELEASE).c_str());
    QDesktopServices::openUrl(QUrl(documentation_url));
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/widget/sketcher_top_bar.moc"
