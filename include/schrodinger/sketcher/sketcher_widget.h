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
class ModifyAtomsMenu;
class ModifyBondsMenu;
class MolModel;
class NonMolecularObject;
class Scene;
class SelectionContextMenu;
class SketcherModel;
enum class ModelKey;
enum class SceneSubset;

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
     * Paste clipboard content into the scene
     */
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
     * Updates the watermark on user drawing atoms or deleting all
     * atoms from the scene
     */
    void updateWatermark();

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
    BackgroundContextMenu* m_background_context_menu = nullptr;

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
    void connectContextMenu(const ModifyAtomsMenu& menu);
    void connectContextMenu(const ModifyBondsMenu& menu);
    void connectContextMenu(const SelectionContextMenu& menu);
    void connectContextMenu(const BackgroundContextMenu& menu);

    /**
     * Show the relevant context menu given some combination of atoms, bonds,
     * and sgroups as emitted by the scene.
     */
    void showContextMenu(
        QGraphicsSceneMouseEvent* event,
        const std::unordered_set<const RDKit::Atom*>& atoms,
        const std::unordered_set<const RDKit::Bond*>& bonds,
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups);

    /**
     * Override QWidget methods to handle keystrokes
     */
    void keyPressEvent(QKeyEvent* event) override;

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
     * Apply the specified SketcherModel operation to the given targets (i.e.
     * atoms, bonds, etc)
     *
     * @param key The SketcherModel key that corresponds to the action to take
     * @param value The SketcherModel value for key
     * @param atoms The atoms, if any, to apply the action to
     * @param bonds The bonds, if any, to apply the action to
     * @param sgroups The substance groups, if any, to apply the action to
     * @param nonmolecular_objects The non-molecular objects (e.g. a reaction
     * arrow or plus sign), if any, to apply the action to
     */
    void applyModelValuePingToTargets(
        const ModelKey key, const QVariant value,
        const std::unordered_set<const RDKit::Atom*> atoms,
        const std::unordered_set<const RDKit::Bond*> bonds,
        const std::unordered_set<const RDKit::SubstanceGroup*> sgroups,
        const std::unordered_set<const NonMolecularObject*>
            nonmolecular_objects);
};

} // namespace sketcher
} // namespace schrodinger
