#include "schrodinger/sketcher/menu/background_context_menu.h"

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

BackgroundContextMenu::BackgroundContextMenu(SketcherModel* model,
                                             QWidget* parent) :
    AbstractContextMenu(parent),
    m_sketcher_model(model)
{
    // Add actions to the menu
    m_save_image_act = addAction("Save Image...", this,
                                 &BackgroundContextMenu::saveImageRequested);
    m_export_to_file_act =
        addAction("Export to File...", this,
                  &BackgroundContextMenu::exportToFileRequested);
    addSeparator();
    m_flip_horizontal_act =
        addAction("Flip All Horizontal", this,
                  &BackgroundContextMenu::flipHorizontalRequested);
    m_flip_vertical_act =
        addAction("Flip All Vertical", this,
                  &BackgroundContextMenu::flipVerticalRequested);
    addSeparator();
    m_undo_act = addAction("Undo", this, &BackgroundContextMenu::undoRequested);
    m_redo_act = addAction("Redo", this, &BackgroundContextMenu::redoRequested);
    addSeparator();
    m_select_all_act = addAction("Select All", this,
                                 &BackgroundContextMenu::selectAllRequested);

    m_cut_copy_manager = new CutCopyActionManager(this);
    m_cut_copy_manager->setModel(model);
    m_cut_copy_manager->setAlwaysCopyAll(true);
    addAction(m_cut_copy_manager->m_copy_action);
    addMenu(m_cut_copy_manager->m_copy_as_menu);
    addAction("Paste", this, [this]() { pasteRequested(pos()); });

    addSeparator();
    addAction("Clear Sketcher", this, &BackgroundContextMenu::clearRequested);
    connect(m_cut_copy_manager, &CutCopyActionManager::copyRequested, this,
            &BackgroundContextMenu::copyRequested);
    connect(m_cut_copy_manager, &CutCopyActionManager::copyAsImageRequested,
            this, &BackgroundContextMenu::copyAsImageRequested);
}

void BackgroundContextMenu::updateActions()
{
    bool has_contents = !m_sketcher_model->sceneIsEmpty();
    auto [can_undo, can_redo] = m_sketcher_model->getUndoStackData();

    // Enable or disable actions
    m_save_image_act->setEnabled(has_contents);
    m_export_to_file_act->setEnabled(has_contents);
    m_flip_horizontal_act->setEnabled(has_contents);
    m_flip_vertical_act->setEnabled(has_contents);
    m_undo_act->setEnabled(can_undo);
    m_redo_act->setEnabled(can_redo);
    m_select_all_act->setEnabled(has_contents);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/background_context_menu.moc"
