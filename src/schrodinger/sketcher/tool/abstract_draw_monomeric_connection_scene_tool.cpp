#include "schrodinger/sketcher/tool/abstract_draw_monomeric_connection_scene_tool.h"

#include <QGraphicsLineItem>

#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"

namespace schrodinger
{
namespace sketcher
{

AbstractDrawMonomericConnectionSceneTool::
    AbstractDrawMonomericConnectionSceneTool(const Fonts& fonts, Scene* scene,
                                             MolModel* mol_model) :
    // the residue name and chain type are irrelevant since this class will
    // never create a new residue, so just pass dummy values
    AbstractDrawMonomerOrMonomericConnectionSceneTool(
        "", rdkit_extensions::ChainType::CHEM, fonts, scene, mol_model)
{
}

void AbstractDrawMonomericConnectionSceneTool::
    updateColorsAfterBackgroundColorChange(bool is_dark_mode)
{
    AbstractDrawMonomerOrMonomericConnectionSceneTool::
        updateColorsAfterBackgroundColorChange(is_dark_mode);
    m_invalid_drag_pen.setColor(is_dark_mode
                                    ? CONNECTION_TOOL_INVALID_DRAG_COLOR_DARK_BG
                                    : CONNECTION_TOOL_INVALID_DRAG_COLOR);
    m_unbound_ap_label_hover_color =
        is_dark_mode ? CONNECTION_TOOL_AP_LABEL_HOVER_COLOR_DARK_BG
                     : CONNECTION_TOOL_AP_LABEL_HOVER_COLOR;
}

bool AbstractDrawMonomericConnectionSceneTool::clickShouldMutate(
    const RDKit::Atom* /* monomer */,
    const MonomerType /* monomer_type */) const
{
    return false;
}

std::optional<HintFragmentMonomerInfo>
AbstractDrawMonomericConnectionSceneTool::
    getHintFragmentMonomerInfoForDragStart()
{
    if (m_drag_start_monomer_item == nullptr) {
        // the drag started over empty space, so we do nothing (unlike
        // DrawMonomerSceneTool, which would add a new monomer here)
        return std::nullopt;
    } else if (!m_drag_start_ap_model_name.empty()) {
        // the drag started over a monomer and it has an available attachment
        // point
        return createHintFragmentMonomerInfoForHintToOrFromExistingMonomer(
            m_drag_start_monomer_item, m_drag_start_ap_model_name);
    }
    // the drag started over a monomer, but that monomer has no available
    // unbound attachment points so we can't drag from it
    return std::nullopt;
}

DragEndInfo
AbstractDrawMonomericConnectionSceneTool::getDragEndInfoForNonMonomerPos(
    const QPointF& scene_pos) const
{
    return scene_pos;
}

void AbstractDrawMonomericConnectionSceneTool::clearHintFragmentItem()
{
    AbstractDrawMonomerOrMonomericConnectionSceneTool::clearHintFragmentItem();
    delete m_invalid_drag_item;
    m_invalid_drag_item = nullptr;
}

void AbstractDrawMonomericConnectionSceneTool::createHintFragmentItem(
    HintFragmentMonomerInfo& monomer_one_info,
    HintFragmentMonomerInfo& monomer_two_info)
{
    if (monomer_one_info.atom_idx != NEW_MONOMER_FROM_INVALID_DRAG &&
        monomer_two_info.atom_idx != NEW_MONOMER_FROM_INVALID_DRAG) {
        AbstractDrawMonomerOrMonomericConnectionSceneTool::
            createHintFragmentItem(monomer_one_info, monomer_two_info);
        return;
    }
    auto start_qcoords = to_scene_xy(monomer_one_info.pos);
    auto end_qcoords = to_scene_xy(monomer_two_info.pos);
    auto* line_item = new QGraphicsLineItem({start_qcoords, end_qcoords});
    m_invalid_drag_item = line_item;
    line_item->setPen(m_invalid_drag_pen);
    m_scene->addItem(line_item);
    line_item->setVisible(true);
}

void AbstractDrawMonomericConnectionSceneTool::addDragStructureToMolModel(
    const HintFragmentMonomerInfo& hint_start_monomer_info,
    const HintFragmentMonomerInfo& hint_end_monomer_info)
{
    if (hint_start_monomer_info.atom_idx != NEW_MONOMER_FROM_INVALID_DRAG &&
        hint_end_monomer_info.atom_idx != NEW_MONOMER_FROM_INVALID_DRAG) {
        AbstractDrawMonomerOrMonomericConnectionSceneTool::
            addDragStructureToMolModel(hint_start_monomer_info,
                                       hint_end_monomer_info);
    }
}

} // namespace sketcher
} // namespace schrodinger
