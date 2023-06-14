#include "schrodinger/sketcher/tool/draw_r_group_scene_tool.h"

#include <algorithm>
#include <unordered_set>

#include <GraphMol/ROMol.h>

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"

namespace schrodinger
{
namespace sketcher
{

AbstractDrawRGroupSceneTool::AbstractDrawRGroupSceneTool(Scene* scene,
                                                         MolModel* mol_model) :
    AbstractDrawAtomSceneTool(scene, mol_model)
{
}

void AbstractDrawRGroupSceneTool::mutateAtom(const RDKit::Atom* const atom)
{
    m_mol_model->mutateRGroup(atom, getRGroupNums(1)[0]);
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

DrawRGroupSceneTool::DrawRGroupSceneTool(const int r_group_num, Scene* scene,
                                         MolModel* mol_model) :
    AbstractDrawRGroupSceneTool(scene, mol_model),
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
    return get_r_group_number(atom) == m_r_group_num;
}

DrawIncrementingRGroupSceneTool::DrawIncrementingRGroupSceneTool(
    Scene* scene, MolModel* mol_model) :
    AbstractDrawRGroupSceneTool(scene, mol_model)
{
}

std::vector<unsigned int>
DrawIncrementingRGroupSceneTool::getRGroupNums(const size_t how_many) const
{
    return m_mol_model->getNextRGroupNumbers(how_many);
}

bool DrawIncrementingRGroupSceneTool::shouldDrawBondForClickOnAtom(
    const RDKit::Atom* const atom) const
{
    return false;
}

} // namespace sketcher
} // namespace schrodinger
