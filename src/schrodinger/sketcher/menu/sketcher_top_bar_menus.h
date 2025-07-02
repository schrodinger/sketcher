#pragma once
#include <QAction>
#include <QMenu>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{
class CutCopyActionManager;
class SketcherModel;

class SKETCHER_API ImportMenu : public QMenu
{
    Q_OBJECT
  public:
    ImportMenu(QWidget* parent = nullptr);
    QAction* m_import_from_file_act = nullptr;
    QAction* m_paste_in_text_act = nullptr;
    QAction* m_replace_content_act = nullptr;
};

class SKETCHER_API ExportMenu : public QMenu
{
    Q_OBJECT
  public:
    ExportMenu(QWidget* parent = nullptr);
    QAction* m_save_image_act = nullptr;
    QAction* m_export_to_file_act = nullptr;
};

class SKETCHER_API MoreActionsMenu : public QMenu
{
    Q_OBJECT
  public:
    MoreActionsMenu(SketcherModel* model, QWidget* parent = nullptr);
    void setModel(SketcherModel* model);

    QAction* m_undo_act = nullptr;
    QAction* m_redo_act = nullptr;
    QAction* m_select_all_act = nullptr;
    QAction* m_clear_selection_act = nullptr;
    QMenu* m_expand_selection_menu = nullptr;
    QAction* m_invert_selection_act = nullptr;
    CutCopyActionManager* m_cut_copy_manager = nullptr;
    QAction* m_paste_act = nullptr;
    QMenu* m_modify_structure_menu = nullptr;
    QAction* m_add_custom_fragment_act = nullptr;
    QAction* m_fit_to_screen_act = nullptr;

    QAction* m_flip_horizontal_act = nullptr;
    QAction* m_flip_vertical_act = nullptr;
    QAction* m_aromatize_act = nullptr;
    QAction* m_kekulize_act = nullptr;
    QAction* m_add_explicit_hydrogens_act = nullptr;
    QAction* m_remove_explicit_hydrogens_act = nullptr;
};

class SKETCHER_API ConfigureViewMenu : public QMenu
{
    Q_OBJECT
  public:
    ConfigureViewMenu(QWidget* parent = nullptr);
    QAction* m_valence_errors_act = nullptr;
    QAction* m_heteroatom_colors_act = nullptr;
    QAction* m_stereo_labels_act = nullptr;
    QAction* m_implicit_hydrogens_act = nullptr;
    QAction* m_preferences_act = nullptr;
};

class SKETCHER_API HelpMenu : public QMenu
{
    Q_OBJECT
  public:
    HelpMenu(QWidget* parent = nullptr);
    QAction* m_help_act = nullptr;
  protected slots:
    void showWelcomeDialog();
    void showAbout2DSketcherDialog();
};

} // namespace sketcher
} // namespace schrodinger
