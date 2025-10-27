#pragma once

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/menu/abstract_context_menu.h"

namespace schrodinger
{

namespace rdkit_extensions
{
enum class Format;
}
namespace sketcher
{

class CutCopyActionManager;
class SketcherModel;
enum class SceneSubset;

class SKETCHER_API BackgroundContextMenu : public AbstractContextMenu
{
    Q_OBJECT
  public:
    BackgroundContextMenu(SketcherModel* model, QWidget* parent = nullptr);

  signals:
    void saveImageRequested();
    void exportToFileRequested();
    void undoRequested();
    void redoRequested();
    void flipHorizontalRequested();
    void flipVerticalRequested();
    void selectAllRequested();
    void copyRequested(rdkit_extensions::Format format, SceneSubset subset);
    void copyAsImageRequested();
    void pasteRequested(QPoint pos);
    void clearRequested();

  protected:
    SketcherModel* m_sketcher_model = nullptr;

    QAction* m_save_image_act = nullptr;
    QAction* m_export_to_file_act = nullptr;
    QAction* m_flip_horizontal_act = nullptr;
    QAction* m_flip_vertical_act = nullptr;
    QAction* m_undo_act = nullptr;
    QAction* m_redo_act = nullptr;
    QAction* m_select_all_act = nullptr;
    CutCopyActionManager* m_cut_copy_manager = nullptr;

    /**
     * Update actions to be enabled or checked depending on the state of the
     * model.
     */
    virtual void updateActions();
};

} // namespace sketcher
} // namespace schrodinger
