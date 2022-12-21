#pragma once

#include <memory>
#include <string>
#include <QWidget>
#include "schrodinger/sketcher/definitions.h"

namespace Ui
{
class SketcherBottomBarForm;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Bottom bar for the sketcher widget meant for use inside Maestro.
 */
class SKETCHER_API SketcherBottomBar : public QWidget
{
    Q_OBJECT

  public:
    SketcherBottomBar(QWidget* parent = nullptr);
    ~SketcherBottomBar();

  signals:
    void resetRequested();
    void reloadRequested();
    void saveAsNewEntryRequested(std::string entry_title);
    void saveChangesRequested();

  private:
    std::unique_ptr<Ui::SketcherBottomBarForm> ui;

  private slots:
    /**
     * Respond to click on "Save as new..." by showing a file dialog.
     *
     * If the user accepts a file, emit `saveAsNewEntryRequested` with the
     * selected file path as its parameter.
     */
    void onSaveAsNewClicked();
};

} // namespace sketcher
} // namespace schrodinger
