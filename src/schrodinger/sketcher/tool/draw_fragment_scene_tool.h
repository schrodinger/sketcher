#pragma once

#include <memory>
#include <string>

#include <GraphMol/ROMol.h>

#include <QList>
#include <QGraphicsItem>
#include <QGraphicsItemGroup>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/atom_item_settings.h"
#include "schrodinger/sketcher/molviewer/bond_item_settings.h"
#include "schrodinger/sketcher/tool/abstract_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{

class Fonts;

/**
 * A graphics item that follows the mouse cursor and shows where the fragment
 * would be added if the user were to click.
 */
class HintFragmentItem : public QGraphicsItemGroup
{
  public:
    /**
     * @param fragment The fragment to display.  Note that the attachment point
     * will not be painted.
     * @param fonts The fonts to use for displaying the fragment.  This object
     * must not be destroyed while this graphics item is in use.
     * @param atom_item_settings The settings for displaying atom items.  This
     * object will be copied.  Note that the color scheme will be ignored, as
     * hints are always displayed in blue.
     * @param bond_item_settings The settings for displaying bond items.  This
     * object will be copied.  Note that the bond width and color will be
     * ignored, as hints are always displayed in blue.
     * @param parent The parent graphics item, if any.
     */
    HintFragmentItem(const RDKit::ROMol& fragment, const Fonts& fonts,
                     const AtomItemSettings& atom_item_settings,
                     const BondItemSettings& bond_item_settings,
                     QGraphicsItem* parent = nullptr);
    void updateConformer(const RDKit::Conformer& conformer);

  protected:
    RDKit::ROMol m_frag;
    /// A list of all child AtomItems
    QList<QGraphicsItem*> m_atom_items;
    /// A list of all child BondItems
    QList<QGraphicsItem*> m_bond_items;
    AtomItemSettings m_atom_item_settings;
    BondItemSettings m_bond_item_settings;
};

/**
 * A scene tool for drawing molecule fragments.  These fragments must have
 * exactly one attachment point, showing where the fragment will be attached.
 */
class SKETCHER_API DrawFragmentSceneTool : public AbstractSceneTool
{
  public:
    /**
     * @param fragment The fragment to draw.  This fragments must have exactly
     * one attachment point
     * @param fonts The fonts to use for displaying the fragment hint.  This
     * object must not be destroyed while this scene tool is in use.
     * @param atom_item_settings The settings for displaying atom item hints.
     * This object will be copied.  Note that the color scheme will be ignored,
     * as hints are always displayed in blue.
     * @param bond_item_settings The settings for displaying bond item hints.
     * This object will be copied.  Note that the bond width and color will be
     * ignored, as hints are always displayed in blue.
     * @param scene The scene that this tool will be used with
     * @param mol_model The MolModel used in scene
     */
    DrawFragmentSceneTool(const RDKit::ROMol& fragment, const Fonts& fonts,
                          const AtomItemSettings& atom_item_settings,
                          const BondItemSettings& bond_item_settings,
                          Scene* scene, MolModel* mol_model);

    // Overridden AbstractSceneTool methods
    std::vector<QGraphicsItem*> getGraphicsItems() override;
    void onMouseMove(QGraphicsSceneMouseEvent* const event) override;
    void onMouseLeave() override;

  protected:
    RDKit::ROMol m_frag;
    HintFragmentItem m_hint_item;
};

/**
 * @return a scene tool for drawing the specified ring.  See the
 * DrawFragmentSceneTool constructor for parameter documentation.
 */
std::shared_ptr<DrawFragmentSceneTool>
get_draw_fragment_scene_tool(const RingTool& ring_tool, const Fonts& fonts,
                             const AtomItemSettings& atom_item_settings,
                             const BondItemSettings& bond_item_settings,
                             Scene* scene, MolModel* mol_model);

/**
 * @return a scene tool for drawing the specified molecule.  See the
 * DrawFragmentSceneTool constructor for parameter documentation.
 *
 * @overload
 */
std::shared_ptr<DrawFragmentSceneTool>
get_draw_fragment_scene_tool(const std::string& text_mol, const Fonts& fonts,
                             const AtomItemSettings& atom_item_settings,
                             const BondItemSettings& bond_item_settings,
                             Scene* scene, MolModel* mol_model);

/**
 * @return the appropriate SMILES string for the specified ring tool
 */
std::string get_smiles_for_ring_tool(const RingTool& ring_tool);

} // namespace sketcher
} // namespace schrodinger
