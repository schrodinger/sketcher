#pragma once

#include <memory>

#include <boost/shared_ptr.hpp>
#include <QWidget>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/rdkit_extensions/convert.h"

class QGraphicsPixmapItem;
class QUndoStack;

namespace RDKit
{
class ROMol;
class ChemicalReaction;
} // namespace RDKit

namespace Ui
{
class SketcherWidgetForm;
}

namespace schrodinger
{

namespace sketcher
{

class MolModel;
class Scene;
class SketcherModel;
enum class ModelKey;

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
};

} // namespace sketcher
} // namespace schrodinger
