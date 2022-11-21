/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <string>

#include <QtGlobal>
#include <QGraphicsScene>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/molviewer/fonts.h"

class QObject;
class QFont;

namespace RDKit
{
class ROMol;
}

namespace schrodinger
{

namespace rdkit_extensions
{
enum class Format;
}

namespace sketcher
{

/**
 * A Qt graphics scene for displaying molecules.
 */
class SKETCHER_API Scene : public QGraphicsScene
{
    Q_OBJECT
  public:
    Scene(QObject* parent = nullptr);

    // TODO: mimic sketcherScene interface for importing/exporting

    /**
     * Load an RDKit molecule into this scene
     *
     * @param mol The RDKit molecule to load.  Note that this class will store a
     * copy of this molecule and will not take ownership of the passed-in
     * instance.
     */
    void loadMol(const RDKit::ROMol& mol);
    /**
     * Load an RDKit molecule into this scene
     *
     * @param mol A shared pointer to the RDKit molecule to load.  After this
     * method is called, the molecule (including the conformation) must no
     * longer be modified by the calling code.
     */
    void loadMol(std::shared_ptr<RDKit::ROMol> mol);

    /**
     * Create a molecule from a text string and load that into the scene.
     * Atomic coordinates will be automatically generated using coordgen.
     *
     * @param text input data to load
     * @param format format to parse
     */
    void importText(const std::string& text,
                    schrodinger::rdkit_extensions::Format format);

    std::shared_ptr<RDKit::ROMol> getRDKitMolecule() const;

    // Getters and setters for changing settings
    qreal fontSize() const;
    void setFontSize(qreal size);

    bool allAtomsShown() const;
    void setAllAtomsShown(bool value);

    CarbonLabels carbonsLabeled() const;
    void setCarbonsLabeled(CarbonLabels value);

    bool valenceErrorsShown() const;
    void setValenceErrorsShown(bool value);

    qreal bondWidth() const;
    void setBondWidth(qreal value);

    qreal doubleBondSpacing() const;
    void setDoubleBondSpacing(qreal value);

  protected:
    std::shared_ptr<RDKit::ROMol> m_mol;
    Fonts m_fonts;
    AtomItemSettings m_atom_item_settings;
    BondItemSettings m_bond_item_settings;

    /**
     * Call updateCachedData() on all AtomItems and BondItems in the scene.
     * (BondItems always need updating after their bound AtomItems are modified
     * in any way.)
     */
    void updateAtomAndBondItems();

    /**
     * Call updateCachedData() on all BondItems (but not AtomItems) in the
     * scene.
     */
    void updateBondItems();

  private:
    /**
     * @return window for the parent widget
     */
    QWidget* window() const;

    /**
     * Overrides the drag/drop events to import text files of supported formats
     */
    void dragEnterEvent(QGraphicsSceneDragDropEvent* event) override;
    void dragMoveEvent(QGraphicsSceneDragDropEvent* event) override;
    void dragLeaveEvent(QGraphicsSceneDragDropEvent* event) override;
    void dropEvent(QGraphicsSceneDragDropEvent* event) override;
};

} // namespace sketcher
} // namespace schrodinger
