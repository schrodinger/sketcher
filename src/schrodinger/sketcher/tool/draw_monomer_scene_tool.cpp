#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"

#include <cmath>

#include <QtMath>
#include <QPen>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/scene.h"

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
    m_fonts(&fonts)
{
    m_highlight_types = InteractiveItemFlag::MONOMER;
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
    // TODO: use a picture of the monomer with the outline so that the cursor
    //       hint indicates the monomer type
    return render_text_to_pixmap(QString::fromStdString(m_res_name),
                                 m_fonts->m_cursor_hint_font,
                                 CURSOR_HINT_COLOR);
}

} // namespace sketcher
} // namespace schrodinger
