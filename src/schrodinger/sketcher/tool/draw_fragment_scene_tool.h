#pragma once

#include <memory>
#include <optional>
#include <string>
#include <utility>

#include <rdkit/GraphMol/ROMol.h>

#include <QList>
#include <QGraphicsItem>
#include <QGraphicsItemGroup>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/atom_display_settings.h"
#include "schrodinger/sketcher/molviewer/bond_display_settings.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

class QPointF;

namespace schrodinger
{
namespace sketcher
{

class Fonts;

/**
 * A graphics item for displaying the fragment
 */
class AbstractHintFragmentItem : public QGraphicsItemGroup
{
  public:
    /**
     * @param fragment The fragment to display.  Note that the attachment point
     * will not be painted.
     * @param fonts The fonts to use for displaying the fragment.  This object
     * must not be destroyed while this graphics item is in use.
     * @param atom_display_settings The settings for displaying atom items. This
     * object will be copied.  Note that the color scheme will be ignored, as
     * hints are always displayed in blue.
     * @param bond_display_settings The settings for displaying bond items. This
     * object will be copied.  Note that the bond width and color will be
     * ignored, as hints are always displayed in blue.
     * @param parent The parent graphics item, if any.
     */
    AbstractHintFragmentItem(const RDKit::ROMol& fragment, const Fonts& fonts,
                             const AtomDisplaySettings& atom_display_settings,
                             const BondDisplaySettings& bond_display_settings,
                             QGraphicsItem* parent = nullptr);

    /**
     * Create a copy of the specified graphics item.  Note that the copy will
     * *not* be parented and will not be added to any scene.
     */
    AbstractHintFragmentItem(const AbstractHintFragmentItem& frag_item);

  protected:
    RDKit::ROMol m_frag;
    const Fonts* m_fonts = nullptr;
    /// A list of all child AtomItems
    QList<QGraphicsItem*> m_atom_items;
    /// A list of all child BondItems
    QList<QGraphicsItem*> m_bond_items;
    /// A list of all child SGroupItems
    QList<QGraphicsItem*> m_s_group_items;
    AtomDisplaySettings m_atom_display_settings;
    BondDisplaySettings m_bond_display_settings;

    /**
     * Create all graphics items required to represent the fragment
     */
    void createGraphicsItems();

    /**
     * Update the atom item and bond item settings (e.g. color scheme, bond line
     * width) to override the passed-in settings.  Subclasses may override this
     * method to change additional settings.
     */
    virtual void updateSettings();
};

/**
 * A graphics item that that shows the fragment in blue at same scale as the
 * core.  This graphics item is used to, e.g., display where the fragment would
 * be added if the user were to click on an atom or bond.
 */
class HintFragmentItem : public AbstractHintFragmentItem
{
  public:
    HintFragmentItem(const RDKit::ROMol& fragment, const Fonts& fonts,
                     const AtomDisplaySettings& atom_display_settings,
                     const BondDisplaySettings& bond_display_settings,
                     QGraphicsItem* parent = nullptr);
    void updateConformer(const RDKit::Conformer& conformer);

    /**
     * Convert all specified bonds to single bonds, and convert all other bonds
     * back to their original type (in case they were changed by a previous call
     * to this method).  This is used to avoid valence errors where the fragment
     * is joined to the core.
     */
    void updateSingleBondMutations(
        const std::unordered_set<RDKit::Bond*>& bonds_to_mutate);

  protected:
    std::vector<RDKit::Bond::Bond::BondType> m_orig_bond_types;
    void updateSettings() override;
};

/**
 * A graphics item that renders the fragment using extra wide lines, suitable
 * for rendering the cursor hint pixmap.  (The wider lines are required because
 * the pixmap is only 28x28 pixels.)
 */
class HintPixmapFragmentItem : public AbstractHintFragmentItem
{
  public:
    HintPixmapFragmentItem(const RDKit::ROMol& fragment, const Fonts& fonts,
                           const AtomDisplaySettings& atom_display_settings,
                           const BondDisplaySettings& bond_display_settings,
                           QGraphicsItem* parent = nullptr);
    HintPixmapFragmentItem(const AbstractHintFragmentItem& frag_item);

  protected:
    void updateSettings() override;
};

/**
 * A scene tool for drawing molecule fragments.  These fragments must have
 * exactly one attachment point, showing where the fragment will be attached.
 */
class SKETCHER_API DrawFragmentSceneTool : public StandardSceneToolBase
{
  public:
    /**
     * @param fragment The fragment to draw.  This fragments must have exactly
     * one attachment point
     * @param fonts The fonts to use for displaying the fragment hint.  This
     * object must not be destroyed while this scene tool is in use.
     * @param atom_display_settings The settings for displaying atom item hints.
     * This object will be copied.  Note that the color scheme will be ignored,
     * as hints are always displayed in blue.
     * @param bond_display_settings The settings for displaying bond item hints.
     * This object will be copied.  Note that the bond width and color will be
     * ignored, as hints are always displayed in blue.
     * @param scene The scene that this tool will be used with
     * @param mol_model The MolModel used in scene
     */
    DrawFragmentSceneTool(const RDKit::ROMol& fragment, const Fonts& fonts,
                          const AtomDisplaySettings& atom_display_settings,
                          const BondDisplaySettings& bond_display_settings,
                          Scene* scene, MolModel* mol_model);

    // Overridden AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    QPixmap createDefaultCursorPixmap() const override;
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void onMouseLeave() override;
    void onLeftButtonPress(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonRelease(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonClick(QGraphicsSceneMouseEvent* const event) override;
    void onLeftButtonDragMove(QGraphicsSceneMouseEvent* const event) override;
    void
    onLeftButtonDragRelease(QGraphicsSceneMouseEvent* const event) override;

  protected:
    RDKit::ROMol m_frag;
    HintFragmentItem m_hint_item;
    bool m_cursor_pixmap_shown = true;

    /**
     * Determine the appropriate fragment conformer for a given set of Scene
     * coordinates.
     * @param scene_pos The Scene coordinates
     * @return A tuple of
     *   - The fragment conformer
     *   - The existing atom (if any) that has the same coordinates as the
     *     fragment's attachment point parent atom (i.e. the atom bound to the
     *     attachment point dummy atom).  Nullptr otherwise.
     */
    std::pair<RDKit::Conformer, const RDKit::Atom*>
    getFragConfAndCoreAtomForScenePos(const QPointF& scene_pos) const;

    /**
     * Determine the appropriate fragment conformer for a click-and-drag to the
     * given set of Scene coordinates.  The start of the drag is taken from
     * m_mouse_press_scene_pos.
     * @param scene_pos The Scene coordinates
     * @return The fragment conformer if the drag was started over empty space.
     * nullopt otherwise.
     */
    std::optional<RDKit::Conformer>
    getConformerForDragToScenePos(const QPointF& scene_pos) const;

    /**
     * Add the fragment to MolModel using the given conformer.  Also hide the
     * fragment hint, as it will overlay the newly added fragment itself.
     * @param conf The fragment conformer
     * @param overlay_atom The atom that the mouse cursor was over (or bond
     * starting atom if the cursor was over a bond).  This will be used to
     * determine which fragment atoms should replace core atoms.  Should be
     * nullptr if the cursor was not over the molecule.
     */
    void addFragToModel(const RDKit::Conformer& conf,
                        const RDKit::Atom* const overlay_atom = nullptr);
};

/**
 * @return a scene tool for drawing the specified ring.  See the
 * DrawFragmentSceneTool constructor for additional parameter documentation.
 *
 * @param ring_tool The ring to use for the fragment
 */
std::shared_ptr<DrawFragmentSceneTool>
get_draw_fragment_scene_tool(const RingTool& ring_tool, const Fonts& fonts,
                             const AtomDisplaySettings& atom_display_settings,
                             const BondDisplaySettings& bond_display_settings,
                             Scene* scene, MolModel* mol_model);

/**
 * @return a scene tool for drawing the specified molecule.  See the
 * DrawFragmentSceneTool constructor for additional parameter documentation.
 *
 * @param text_mol The molecule to use for the fragment
 * @param format The format that text_mol is in
 *
 * @overload
 */
std::shared_ptr<DrawFragmentSceneTool>
get_draw_fragment_scene_tool(const std::string& text_mol, const Fonts& fonts,
                             const AtomDisplaySettings& atom_display_settings,
                             const BondDisplaySettings& bond_display_settings,
                             Scene* scene, MolModel* mol_model,
                             const rdkit_extensions::Format format =
                                 rdkit_extensions::Format::AUTO_DETECT);

/**
 * @return the appropriate SMILES string for the specified ring tool
 */
std::string get_smiles_for_ring_tool(const RingTool& ring_tool);

} // namespace sketcher
} // namespace schrodinger
