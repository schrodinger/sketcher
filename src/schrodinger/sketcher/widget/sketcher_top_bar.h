#pragma once

#include <memory>
#include <string>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/widget/sketcher_view.h"

class QWidget;

namespace Ui
{
class SketcherTopBar;
}

namespace schrodinger
{
namespace rdkit_extensions
{
enum class Format;
}
namespace sketcher
{
class ConfigureViewMenu;
class ExportMenu;
class HelpMenu;
class ImportMenu;
class MoreActionsMenu;
enum class SceneSubset;

/**
 * Top bar for the sketcher widget.
 */
class SKETCHER_API SketcherTopBar : public SketcherView
{
    Q_OBJECT

  public:
    SketcherTopBar(QWidget* parent = nullptr);
    ~SketcherTopBar();
    void setModel(SketcherModel* model) override;

    /** Check the given key event and if it matches with any shortcut
     *  defined in the topbar, manually trigger it.
     *  If event is handled by this method, returns true and accept the
     *  event, otherwise return false.
     */
    bool handleShortcutAction(const QKeySequence& key_seq);

    MoreActionsMenu* m_more_actions_menu = nullptr;

  signals:
    void clearSketcherRequested();
    void importTextRequested(const std::string& text,
                             const rdkit_extensions::Format format);
    void saveImageRequested();
    void exportToFileRequested();
    void undoRequested();
    void redoRequested();
    void cleanupRequested(bool selection_only);
    void selectAllRequested();
    void clearSelectionRequested();
    void expandSelectionRequested();
    void invertSelectionRequested();
    void cutRequested(rdkit_extensions::Format format);
    void copyRequested(rdkit_extensions::Format format, SceneSubset subset);
    void pasteRequested();
    void flipHorizontalRequested();
    void flipVerticalRequested();
    void aromatizeRequested();
    void kekulizeRequested();
    void addExplicitHydrogensRequested();
    void removeExplicitHydrogensRequested();
    void addCustomFragmentRequested(const std::string& text);
    void adjustRenderingSettingsRequested();
    void fitToScreenRequested(bool selection_only = false);

  protected:
    std::unique_ptr<Ui::SketcherTopBar> ui;
    ImportMenu* m_import_menu = nullptr;
    ExportMenu* m_export_menu = nullptr;
    ConfigureViewMenu* m_configure_view_menu = nullptr;
    HelpMenu* m_help_menu = nullptr;

    /**
     * Generate a vector that identifies setters associated with each model key.
     *
     * This function is meant to be run a single time during instantiation in
     * order to assign relevant values to `m_setter_packets` and
     * `m_signal_packets`.
     */
    void generatePackets() override;

  protected slots:
    void updateWidgetsEnabled() override;

  private:
    /**
     * Assign button IDs and connect signals.
     */
    void initButtons();

    /**
     * Create menus, assign them to buttons, and connect their signals.
     */
    void initMenus();

  private slots:

    /**
     * Respond to a click on the "Import from File..." action.
     *
     * Present the user with a file selection dialog. If accepted, emit a
     * `importTextRequested` signal with the specified file path.
     */
    void onImportFromFileClicked();

    /**
     * Respond to a click on the "Paste in Text..." action.
     *
     * Present the user with a text box dialog. If accepted, emit a
     * `importTextRequested` signal with the specified text.
     */
    void onPasteInTextClicked();

    /**
     * Respond to a click on the "Add Custom Fragment..." action.
     *
     * Present the user with the custom fragment dialog. If accepted, emit a
     * `addCustomFragmentRequested` signal with information about the specified
     * fragment.
     */
    void onAddCustomFragmentClicked();

    /**
     * Respond to a click on the "Help..." action, opening the documentation
     */
    virtual void onHelpClicked();
};

} // namespace sketcher
} // namespace schrodinger
