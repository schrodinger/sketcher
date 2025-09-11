#pragma once

#include <memory>
#include <unordered_set>

#include <boost/shared_ptr.hpp>
#include <QWidget>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/rdkit_extensions/convert.h"

class QGraphicsPixmapItem;
class QGraphicsSceneMouseEvent;
class QUndoStack;

namespace RDKit
{
class Atom;
class Bond;
class ROMol;
class ChemicalReaction;
class SubstanceGroup;
} // namespace RDKit

namespace Ui
{
class SketcherWidgetForm;
}

namespace schrodinger
{

namespace sketcher
{

class AtomContextMenu;
class BackgroundContextMenu;
class BondContextMenu;
class BracketSubgroupContextMenu;
class ModifyAtomsMenu;
class ModifyBondsMenu;
class MolModel;
class NonMolecularObject;
class Scene;
class SelectionContextMenu;
class SketcherModel;
enum class ImageFormat;
enum class ModelKey;
enum class SceneSubset;
enum class ColorScheme;
class FileExportDialog;
class FileSaveImageDialog;
class RenderingSettingsDialog;
class ModelObjsByType;

/**
 * Sketcher widget meant for use in LiveDesign and Maestro.
 */
class SKETCHER_API SketcherWidget : public QWidget
{
    Q_OBJECT

  public:
    SketcherWidget(QWidget* parent = nullptr);
    ~SketcherWidget();

    /**
     * @param mol/rxn molecule or reaction to add to the the current scene
     */
    void addRDKitMolecule(const RDKit::ROMol& mol);
    void addRDKitReaction(const RDKit::ChemicalReaction& rxn);

    /**
     * @return the sketcher contents as either a RDKit molecule or reaction
     */
    boost::shared_ptr<RDKit::ROMol> getRDKitMolecule() const;
    boost::shared_ptr<RDKit::ChemicalReaction> getRDKitReaction() const;

    /**
     * @param text serialized molecule/reaction to load into the sketcher
     * @param format specified format of the input text
     */
    void addFromString(const std::string& text,
                       rdkit_extensions::Format format =
                           rdkit_extensions::Format::AUTO_DETECT);

    /**
     * @return the sketcher contents in the request serialized format
     */
    std::string getString(rdkit_extensions::Format format) const;

    /**
     * @return the sketcher contents as an image in the requested format
     */
    QByteArray getImageBytes(ImageFormat format) const;

    /**
     * Clear the scene in its entirety
     */
    void clear();

    /**
     * @return true if the scene is empty
     */
    bool isEmpty() const;

    /**
     * Enable or disable select-only mode, which removes the toolbars and limits
     * user interaction to selecting items.
     */
    void setSelectOnlyMode(bool select_only_mode_enabled);

    /**
     * Set the color scheme, which controls the atom and bond colors as well as
     * the background color.
     */
    void setColorScheme(ColorScheme color_scheme);

    /**
     * Pass the specified keypress to the top bar.  Explicit calls to this
     * method are necessary on Mac to prevent the Maestro window from stealing
     * keypresses.
     */
    bool handleShortcutAction(const QKeySequence& key);

  signals:

    /**
     * Emitted when an atom or bond is added to, modified, or removed from the
     * molecule shown in the Sketcher.  This signal is *not* emitted, however,
     * when the two-dimensional coordinates are changed, or when reaction arrows
     * or plus signs are added or removed.
     *
     * In other words, this signal will be emitted when the SMILES string for
     * the molecule changes (but not necessarily emitted when the *reaction*
     * SMILES string changes, since that can be affected by the position of the
     * reaction arrow).
     */
    void moleculeChanged();

    /**
     * Emitted when any of the following happen
     *   - the two-dimensional coordinates of the molecule are changed
     *   - a reaction arrow or a plus sign is added or removed
     *   - the two-dimensional coordinates of the reaction arrow or plus signs
     *     are changed
     *
     * This signal will *not* be emitted when
     *   -  an atom or bond is added to, modified, or removed from the molecule
     *     (Use the moleculeChanged signal to receive notifications for this.)
     *   - selection changes
     *   - the entire molecule is translated or zoomed.  These actions change
     *     the viewport, but have no effect on the coordinates of the molecule
     *     itself.  Note that rotating the entire molecule *does* affect the
     *     coordinates and *will* cause this signal to be emitted
     *   - display options are changed, such as displaying valence errors,
     *     heteroatom colors, or stereocenter labels
     *   - *during* a click-and-drag operation that affects coordinates, such as
     *     rotating the entire molecule, rotating a portion of the molecule, or
     *     translating a portion of the molecule.  Instead, this signal will
     *     be emitted *at the end* of the click-and-drag.
     */
    void representationChanged();

  protected slots:

    /**
     * Copy the contents of the scene to the clipboard and delete the selection
     *
     * @param format The format of the text structure to copy to clipboard
     */
    void cut(rdkit_extensions::Format format);

    /**
     * Copy the contents of the scene to the clipboard.
     *
     * @param format The format of the text structure to copy to clipboard
     * @param subset The portion of the scene structure to copy to clipboard
     */
    void copy(rdkit_extensions::Format format, SceneSubset subset);

    /**
     * Paste clipboard content into the scene.
     * @param position The position to paste the content at
     */
    void pasteAt(std::optional<QPointF> position);

    void paste();

    /**
     * Import the given text into the scene; optionally clear beforehand
     * depending on the state of the model
     */
    void importText(const std::string& text, rdkit_extensions::Format format);

    /**
     * Present the user with an "Export to File" dialog.
     */
    void showFileExportDialog();

    /**
     * Present the user with a "Save Image" dialog.
     */
    void showFileSaveImageDialog();

    /**
     * Show the rendering settings dialog, which allows the user to change
     * rendering options such as atom colors, bond widths, and background color.
     */
    void showRenderingSettingsDialog();

    /**
     * Show the Edit Atom Properties dialog for the specified atom
     * @param atom The atom to load into the dialog
     * @param set_to_allowed_list Whether to switch the dialog to display an
     * allowed list query
     */
    void showEditAtomPropertiesDialog(const RDKit::Atom* const atom,
                                      const bool set_to_allowed_list);

    /**
     * Updates the watermark on user drawing atoms or deleting all
     * atoms from the scene
     */
    void updateWatermark();

    /**
     * Show the substance group bracket dialog.  If the dialog is accepted, a
     * new S-group will be created for the specified atoms.
     */
    void showBracketSubgroupDialogForAtoms(
        const std::unordered_set<const RDKit::Atom*>& atoms);

    /**
     * Show the substance group bracket dialog.  If the dialog is accepted, the
     * specified S-group will be modified.
     */
    void showBracketSubgroupDialogForSGroup(
        const RDKit::SubstanceGroup* const s_group);

  protected:
    /**
     * Widget UI form and layout
     */
    std::unique_ptr<Ui::SketcherWidgetForm> m_ui;

    /**
     * Models and scene owned by the widget
     */
    QUndoStack* m_undo_stack = nullptr;
    MolModel* m_mol_model = nullptr;
    SketcherModel* m_sketcher_model = nullptr;
    Scene* m_scene = nullptr;

    /**
     * Context menus
     */
    AtomContextMenu* m_atom_context_menu = nullptr;
    BondContextMenu* m_bond_context_menu = nullptr;
    SelectionContextMenu* m_selection_context_menu = nullptr;
    BracketSubgroupContextMenu* m_sgroup_context_menu = nullptr;
    BackgroundContextMenu* m_background_context_menu = nullptr;

    /**
     * Dialogs
     */
    FileExportDialog* m_file_export_dialog = nullptr;
    FileSaveImageDialog* m_file_save_image_dialog = nullptr;
    RenderingSettingsDialog* m_rendering_settings_dialog = nullptr;

    /**
     * Watermark centered in the Scene; only shown when no atoms are present
     */
    QGraphicsPixmapItem* m_watermark_item = nullptr;

    bool m_select_only_mode_active = false;

    /**
     * Methods for getting and setting the system clipboard. These methods are
     * overriden in the unit test since using the real clipboard can lead to
     * intermittent failures on buildbot.
     */
    virtual std::string getClipboardContents() const;
    virtual void setClipboardContents(std::string text) const;

    /**
     *  Connects slots to the model and various widget tool bars
     */
    void connectTopBarSlots();
    void connectSideBarSlots();
    void connectContextMenu(const ModifyAtomsMenu& menu);
    void connectContextMenu(const ModifyBondsMenu& menu);
    void connectContextMenu(const SelectionContextMenu& menu);
    void connectContextMenu(const BracketSubgroupContextMenu& menu);
    void connectContextMenu(const BackgroundContextMenu& menu);

    /**
     * Show the relevant context menu given some combination of atoms, bonds,
     * and sgroups as emitted by the scene.
     */
    void showContextMenu(
        QGraphicsSceneMouseEvent* event,
        const std::unordered_set<const RDKit::Atom*>& atoms,
        const std::unordered_set<const RDKit::Bond*>& bonds,
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
        const std::unordered_set<const NonMolecularObject*>&
            non_molecular_objects);

    /**
     * Show or hide the toolbars
     */
    void setToolbarsVisible(const bool visible);

    /**
     * Override QWidget methods to handle keystrokes
     */
    void keyPressEvent(QKeyEvent* event) override;

    /**
     * Update the SketcherModel when handling a keyboard shortcut
     *
     * @param has_targets whether this shortcut will operate on a specific
     * subset of model objects
     * @param kv_pair SketcherModel key and value to set if there is currently
     * no workspace selection
     * @param kv_pairs SketcherModel kays and values to ping if has_targets is
     * true
     * @param targets the model objects to apply kV_pairs (and kv_pair if
     * there's a selection) to
     */
    void updateModelForKeyboardShortcut(
        const bool has_targets, const std::pair<ModelKey, QVariant>& kv_pair,
        std::unordered_map<ModelKey, QVariant> kv_pairs,
        const ModelObjsByType& targets);

    /**
     * Process keyboard shortcuts that are common to all input modes (atomistic,
     * amino acid, and nucleic acid)
     *
     * @param event the key press event to be handled
     * @param cursor_pos the current cursor position
     * @param targets the model objects that the key press should operate on
     * @return whether the key press was handled
     */
    bool handleCommonKeyboardShortcuts(QKeyEvent* event,
                                       const QPointF& cursor_pos,
                                       const ModelObjsByType& targets);

    /**
     * Process keyboard shortcuts when the atomistic tools are active
     *
     * @note see handleCommonKeyboardShortcuts for param documentation
     */
    void handleAtomisticKeyboardShortcuts(QKeyEvent* event,
                                          const QPointF& cursor_pos,
                                          const ModelObjsByType& targets);

    /**
     * Respond to the user clicking on a toolbar button when there is a
     * selection present.  Note that this method is called whenever the user
     * clicks on a tool button, but it's a no-op unless a selection is present.
     *
     * @param key The SketcherModel key associated with the button that was
     * clicked
     * @param value The current SketcherModel value for key
     */
    void onModelValuePinged(ModelKey key, QVariant value);

    /**
     * Respond to a value changing in the model
     * @param color The new SketcherModel background color
     */
    void onBackgroundColorChanged(const QColor& color);

    /**
     * Apply the specified SketcherModel operation to the given targets (i.e.
     * atoms, bonds, etc)
     *
     * @param key The SketcherModel key that corresponds to the action to take
     * @param value The SketcherModel value for key
     * @param atoms The atoms, if any, to apply the action to
     * @param bonds The bonds, if any, to apply the action to
     * @param sgroups The substance groups, if any, to apply the action to
     * @param non_molecular_objects The non-molecular objects (e.g. a reaction
     * arrow or plus sign), if any, to apply the action to
     */
    void applyModelValuePingToTargets(
        const ModelKey key, const QVariant value,
        const std::unordered_set<const RDKit::Atom*> atoms,
        const std::unordered_set<const RDKit::Bond*> bonds,
        const std::unordered_set<const RDKit::SubstanceGroup*> sgroups,
        const std::unordered_set<const NonMolecularObject*>
            non_molecular_objects);
};

} // namespace sketcher
} // namespace schrodinger
