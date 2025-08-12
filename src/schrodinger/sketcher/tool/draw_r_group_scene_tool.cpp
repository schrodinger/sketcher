#include "schrodinger/sketcher/tool/draw_r_group_scene_tool.h"

#include <algorithm>
#include <unordered_set>

#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

AbstractDrawRGroupSceneTool::AbstractDrawRGroupSceneTool(const Fonts& fonts,
                                                         Scene* scene,
                                                         MolModel* mol_model) :
    AbstractDrawAtomSceneTool(fonts, scene, mol_model)
{
}

void AbstractDrawRGroupSceneTool::mutateAtom(const RDKit::Atom* const atom)
{
    m_mol_model->mutateRGroups({atom}, getRGroupNums(1)[0]);
}

void AbstractDrawRGroupSceneTool::addAtom(const RDGeom::Point3D& pos,
                                          const RDKit::Atom* const bound_to)
{
    m_mol_model->addRGroupChain(getRGroupNums(1), {pos}, bound_to);
}

void AbstractDrawRGroupSceneTool::addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                                                   const RDGeom::Point3D& pos2)
{
    m_mol_model->addRGroupChain(getRGroupNums(2), {pos1, pos2});
}

DrawRGroupSceneTool::DrawRGroupSceneTool(const int r_group_num,
                                         const Fonts& fonts, Scene* scene,
                                         MolModel* mol_model) :
    AbstractDrawRGroupSceneTool(fonts, scene, mol_model),
    m_r_group_num(r_group_num)
{
}

std::vector<unsigned int>
DrawRGroupSceneTool::getRGroupNums(const size_t how_many) const
{
    return std::vector<unsigned int>(how_many, m_r_group_num);
}

bool DrawRGroupSceneTool::shouldDrawBondForClickOnAtom(
    const RDKit::Atom* const atom) const
{
    auto r_group_num = rdkit_extensions::get_r_group_number(atom);
    return r_group_num.has_value() && r_group_num.value() == m_r_group_num;
}

QPixmap DrawRGroupSceneTool::createDefaultCursorPixmap() const
{
    return render_text_to_pixmap(QString("R%1").arg(m_r_group_num),
                                 m_fonts->m_cursor_hint_font,
                                 CURSOR_HINT_COLOR);
}

DrawIncrementingRGroupSceneTool::DrawIncrementingRGroupSceneTool(
    const Fonts& fonts, Scene* scene, MolModel* mol_model) :
    AbstractDrawRGroupSceneTool(fonts, scene, mol_model)
{
}

std::vector<unsigned int>
DrawIncrementingRGroupSceneTool::getRGroupNums(const size_t how_many) const
{
    return get_next_r_group_numbers(m_mol_model->getMol(), how_many);
}

bool DrawIncrementingRGroupSceneTool::shouldDrawBondForClickOnAtom(
    const RDKit::Atom* const atom) const
{
    return false;
}

QPixmap DrawIncrementingRGroupSceneTool::createDefaultCursorPixmap() const
{
    return render_text_to_pixmap("R+", m_fonts->m_cursor_hint_font,
                                 CURSOR_HINT_COLOR);
}

} // namespace sketcher
} // namespace schrodinger
