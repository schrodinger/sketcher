
#include "schrodinger/sketcher/menu/cut_copy_action_manager.h"

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

// Controls the default cut/copy format when not specified
const Format DEFAULT_FORMAT = Format::MDL_MOLV3000;

CutCopyActionManager::CutCopyActionManager(QWidget* parent) : QWidget(parent)
{
    m_cut_action = new QAction("Cut", this);
    m_copy_action = new QAction("Copy", this);
    m_copy_as_menu = new QMenu("Copy As", this);
}

void CutCopyActionManager::setModel(SketcherModel* model)
{
    if (m_sketcher_model != nullptr) {
        throw std::runtime_error(
            "The model has already been set on this view.");
    }

    m_sketcher_model = model;
    connect(m_sketcher_model, &SketcherModel::interactiveItemsChanged, this,
            &CutCopyActionManager::updateActions);
    connect(m_sketcher_model, &SketcherModel::selectionChanged, this,
            &CutCopyActionManager::updateActions);
    connect(m_sketcher_model, &SketcherModel::undoStackDataChanged, this,
            &CutCopyActionManager::updateActions);

    // Postpone connecting actions until the model is set since copy
    // queries the model to determine which subset value to emit
    connect(m_cut_action, &QAction::triggered, this,
            [this]() { emit cutRequested(DEFAULT_FORMAT); });
    connect(m_copy_action, &QAction::triggered, this,
            [this]() { emit copyRequested(DEFAULT_FORMAT, getSubset()); });
    initCopyAsMenu();

    // Initialize
    updateActions();
}

void CutCopyActionManager::setAlwaysCopyAll(bool always_copy_all)
{
    m_always_copy_all = always_copy_all;
}

SceneSubset CutCopyActionManager::getSubset()
{
    if (m_always_copy_all || !m_sketcher_model->hasActiveSelection()) {
        return SceneSubset::ALL;
    } else {
        return SceneSubset::SELECTION;
    }
}

void CutCopyActionManager::initCopyAsMenu()
{
    auto init_menu = [&](const auto& format_list, bool is_reaction_format) {
        for (const auto& format : format_list) {
            // Clang < 16 does not allow capturing structured bindings, so
            // we'll extract "manually".
            auto& fmt = std::get<0>(format);
            auto& label = std::get<1>(format);
            auto slot = [this, fmt]() { emit copyRequested(fmt, getSubset()); };
            auto action = m_copy_as_menu->addAction(
                QString::fromStdString(label), this, slot);
            // set a flag on the action to determine its visibility later
            action->setData(QVariant(is_reaction_format));
        }
    };

    init_menu(get_standard_export_formats(), false);
    init_menu(get_reaction_export_formats(), true);

    // Add a separator and the option to export as an image
    auto separator = m_copy_as_menu->addSeparator();
    auto action = m_copy_as_menu->addAction(
        "Image", this, [this]() { emit copyAsImageRequested(); });
    // export as image should only be allowed for full scene
    m_hide_for_selections.push_back(separator);
    m_hide_for_selections.push_back(action);
}

void CutCopyActionManager::updateActions()
{
    bool has_contents = !m_sketcher_model->sceneIsEmpty();
    bool has_selection = m_sketcher_model->hasActiveSelection();
    m_cut_action->setEnabled(has_contents && has_selection);
    m_copy_action->setEnabled(has_contents);
    m_copy_as_menu->setEnabled(has_contents);

    bool export_selection = getSubset() == SceneSubset::SELECTION;
    if (export_selection) {
        m_copy_action->setText("Copy");
        m_copy_as_menu->setTitle("Copy As");
    } else {
        m_copy_action->setText("Copy All");
        m_copy_as_menu->setTitle("Copy All As");
    }

    // For reaction formats to show, 1) there must be a reaction present
    // and 2) if a selection is to be exported, everything must be selected
    bool show_reaction = m_sketcher_model->hasReaction();
    if (show_reaction && export_selection) {
        show_reaction = m_sketcher_model->allItemsSelected();
    }

    for (auto act : m_copy_as_menu->actions()) {
        if (m_hide_for_selections.contains(act)) {
            // export as image and its separator, only show for full scene
            act->setVisible(!export_selection);
        } else {
            // structure vs reaction formats
            act->setVisible(act->data().toBool() == show_reaction);
        }
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/menu/cut_copy_action_manager.moc"
