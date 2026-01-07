#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"

#include <cmath>
#include <memory>

#include <QtMath>
#include <QPen>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"

namespace schrodinger
{
namespace sketcher
{

DrawMonomerSceneTool::DrawMonomerSceneTool(
    const std::string& res_name, const rdkit_extensions::ChainType chain_type,
    const Fonts& fonts, Scene* scene, MolModel* mol_model) :
    StandardSceneToolBase(scene, mol_model),
    m_res_name(res_name),
    m_chain_type(chain_type),
    m_fonts(fonts)
{
    m_highlight_types = InteractiveItemFlag::MONOMER;
    // make sure that the cursor hint font is more easily readable at small size
    m_fonts.m_main_label_font.setBold(true);
    m_fonts.updateFontMetrics();
}

std::vector<QGraphicsItem*> DrawMonomerSceneTool::getGraphicsItems()
{
    auto items = StandardSceneToolBase::getGraphicsItems();
    // TODO
    return items;
}

void DrawMonomerSceneTool::onMouseMove(QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onMouseMove(event);
    if (m_mouse_pressed) {
        // drag logic is handled in onDragMove
        return;
    }
    // QPointF scene_pos = event->scenePos();
    // auto* item = m_scene->getTopInteractiveItemAt(scene_pos,
    //                                               InteractiveItemFlag::MONOMER);
    // TODO: add connector nubs to hovered monomers
}

void DrawMonomerSceneTool::onLeftButtonClick(
    QGraphicsSceneMouseEvent* const event)
{
    StandardSceneToolBase::onLeftButtonClick(event);
    QPointF scene_pos = event->scenePos();
    auto* item = m_scene->getTopInteractiveItemAt(scene_pos,
                                                  InteractiveItemFlag::MONOMER);
    if (item == nullptr) {
        // the click was on empty space, so create a new monomer here
        auto mol_pos = to_mol_xy(scene_pos);
        m_mol_model->addMonomer(m_res_name, m_chain_type, mol_pos);
    } else {
        // the click was on an existing monomer
        //  auto* monomer_item = dynamic_cast<AbstractMonomerItem*>(item);
        //  const auto* monomer = monomer_item->getAtom();
        // TODO
    }
}

QPixmap DrawMonomerSceneTool::createDefaultCursorPixmap() const
{
    // the specific number used here (the "1") doesn't matter - we just need any
    // number to form a proper chain ID
    auto chain_id = rdkit_extensions::toString(m_chain_type) + "1";
    auto monomer =
        rdkit_extensions::makeMonomer(m_res_name, chain_id, 1, false);

    std::shared_ptr<AbstractMonomerItem> monomer_item;
    monomer_item.reset(get_monomer_graphics_item(monomer.get(), m_fonts));
    monomer_item->setMonomerColors(Qt::GlobalColor::transparent,
                                   CURSOR_HINT_COLOR, CURSOR_HINT_COLOR);
    // make sure that the cursor hint is at least a little smaller than an
    // actual monomer
    auto min_scene_size =
        CURSOR_HINT_IMAGE_SIZE * MONOMER_CURSOR_HINT_MIN_SCENE_SIZE_SCALE;
    return cursor_hint_from_graphics_item(monomer_item.get(), min_scene_size);
}

} // namespace sketcher
} // namespace schrodinger
