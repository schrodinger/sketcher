#pragma once

#include <memory>

#include <QGraphicsItem>
#include <QGraphicsItemGroup>
#include <QList>

#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/definitions.h"
#include "schrodinger/sketcher/molviewer/constants.h"

namespace RDKit
{
class RWMol;
}

namespace schrodinger::sketcher
{

class Fonts;
class Scene;

/**
 * A graphics item showing a blue hint structure for monomeric models.
 */
class MonomerHintFragmentItem : public QGraphicsItemGroup
{
  public:
    /**
     * @param fragment The fragment to display. Note that the conformer of this
     * molecule may be modified.
     * @param fonts The fonts to use for displaying the fragment. This object
     * must not be destroyed while this graphics item is in use.
     * @param atom_indices_to_hide The graphics item for these atom will be
     * hidden. Normally used to hide atoms that overlap the existing Sketcher
     * structure.
     * @param bond_index_to_label If >= 0, the attachment points for this
     * connector will be labeled.
     * @param monomer_background_color The color to use for the monomer
     * background. This should normally be the same as the Scene background.
     * Using a transparent color is not recommended, as the lines for monomeric
     * connections will be visible behind the monomer outlines and labels.
     * @param parent The parent graphics item, if any.
     */
    MonomerHintFragmentItem(std::shared_ptr<RDKit::RWMol> fragment,
                            const Fonts& fonts,
                            const std::vector<size_t>& atom_indices_to_hide,
                            const int bond_index_to_label,
                            const QColor monomer_background_color,
                            QGraphicsItem* parent = nullptr);

    /**
     * Rotate the fragment to the specified angle relative to its original
     * conformation.
     *
     * @note: This method doesn't rotate the attachment point labels created
     * when bond_index_to_label is >= 0 (since it's currently only in scenarios
     * where no such labels are drawn).
     */
    void setRotation(const double angle_radians,
                     const int monomer_idx_to_rotate_about);

  protected:
    std::shared_ptr<RDKit::RWMol> m_frag;
    const Fonts* m_fonts = nullptr;
    std::vector<size_t> m_atom_indices_to_hide;
    int m_bond_index_to_label = -1;
    QColor m_monomer_background_color;
    QList<QGraphicsItem*> m_atom_items;
    QList<QGraphicsItem*> m_bond_items;
    RDKit::Conformer m_orig_conf;

    /**
     * Create all graphics items required to represent the fragment
     */
    void createGraphicsItems();

    /**
     * Style all of the graphics items used to represent the fragment
     */
    void styleGraphicsItems();
};

} // namespace schrodinger::sketcher