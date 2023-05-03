#pragma once

#include <memory>

#include <QWidget>

#include "schrodinger/sketcher/definitions.h"

class QGraphicsPixmapItem;
class QUndoStack;

namespace Ui
{
class SketcherWidgetForm;
}

namespace schrodinger
{

namespace rdkit_extensions
{
enum class Format;
}

namespace sketcher
{

class MolModel;
class Scene;
class SketcherModel;

/**
 * Sketcher widget meant for use in LiveDesign and Maestro.
 */
class SKETCHER_API SketcherWidget : public QWidget
{
    Q_OBJECT

  public:
    SketcherWidget(QWidget* parent = nullptr);
    ~SketcherWidget();

  protected:
    /**
     * Import the given text into the scene; optionally clear beforehand
     * depending on the state of the model
     *
     * @param text input data to load into the sketcher
     * @param format format to parse
     */
    void importText(const std::string& text, rdkit_extensions::Format format);

    /**
     * Paste clipboard content into the scene
     */
    void paste();

    /**
     * Present the user with an "Export to File" dialog.
     */
    void showFileExportDialog();

    /**
     * Present the user with a "Save Image" dialog.
     */
    void showFileSaveImageDialog();

    /**
     * Updates the watermark on user drawing atoms or deleting all
     * atoms from the scene
     */
    void updateWatermark();

    std::unique_ptr<Ui::SketcherWidgetForm> m_ui;

    /**
     * Models and scene owned by the widget
     */
    QUndoStack* m_undo_stack = nullptr;
    MolModel* m_mol_model = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
    Scene* m_scene = nullptr;

    /**
     * Watermark centered in the Scene; only shown when no atoms are present
     */
    QGraphicsPixmapItem* m_watermark_item = nullptr;

  private:
    /**
     *  Connects slots to the model and various widget tool bars
     */
    void connectTopBarSlots();
    void connectSideBarSlots();
};

} // namespace sketcher
} // namespace schrodinger
