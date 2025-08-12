#pragma once

#include <QDialog>
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/image_constants.h"
#include "schrodinger/sketcher/dialog/modal_dialog.h"

namespace Ui
{
class RenderingSettingsDialog;
}

namespace schrodinger
{
namespace sketcher
{

struct RenderingSettings;

class SKETCHER_API RenderingSettingsDialog : public ModalDialog
{
    Q_OBJECT
  public:
    RenderingSettingsDialog(SketcherModel* model, QWidget* parent);
    ~RenderingSettingsDialog();

  protected slots:
    void onValuesChanged();

    /**
     * @brief Load the default rendering settings
     */
    void loadDefaults();

  private:
    /**
     * @brief Load the current rendering settings from the sketcher model
     */
    void loadSettingsFromModel();

    /**
     * @brief export the current state to the sketcher model
     */
    void exportSettingsToModel();

    /**
     * @brief Update the GUI elements based on the current settings
     */
    void updateWidgets();

    /**
     * @brief Load the color modes combobox options
     */
    void initColorModes();

    /**
     * @brief Get the current state of the rendering settings from the sketcher
     * model
     */
    RenderingSettings
    getSettingsFromModel(const sketcher::SketcherModel* model) const;

    /**
     * @brief Update the model with the given settings
     */
    void loadSettings(RenderingSettings& settings);

    /**
     * @brief Get the current rendering settings
     */
    RenderingSettings getSettingsFromPanel() const;

    std::unique_ptr<Ui::RenderingSettingsDialog> m_ui;
    SketcherModel* m_sketcher_model = nullptr;
    QTimer* m_update_timer = nullptr;
    bool m_freeze_update_from_model = false; // prevent updates from the model
};

} // namespace sketcher
} // namespace schrodinger
