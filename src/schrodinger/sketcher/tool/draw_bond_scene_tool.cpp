#include "schrodinger/sketcher/tool/draw_bond_scene_tool.h"

#include <tuple>
#include <unordered_map>

#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"

namespace schrodinger
{
namespace sketcher
{

AbstractDrawBondSceneTool::AbstractDrawBondSceneTool(Scene* scene,
                                                     MolModel* mol_model) :
    AbstractDrawSceneTool(scene, mol_model)
{
}

bool AbstractDrawBondSceneTool::shouldDrawBondForClickOnAtom(
    const RDKit::Atom* const atom) const
{
    return true;
}

void AbstractDrawBondSceneTool::onBondClicked(const RDKit::Bond* const bond)
{
    if (bondMatches(bond)) {
        if (m_flippable) {
            m_mol_model->flipBond(bond);
        } else if (m_cycle_on_click) {
            cycleBond(bond);
        }
    } else {
        mutateBond(bond);
    }
}

void AbstractDrawBondSceneTool::onEmptySpaceClicked(const RDGeom::Point3D& pos)
{
    RDGeom::Point3D pos2(pos);
    pos2.x += BOND_LENGTH;
    addTwoBoundAtoms(pos, pos2);
}

DrawBondSceneTool::DrawBondSceneTool(BondTool bond_tool, Scene* scene,
                                     MolModel* mol_model) :
    AbstractDrawBondSceneTool(scene, mol_model)
{
    std::tie(m_bond_type, m_bond_dir, m_flippable) =
        BOND_TOOL_BOND_MAP.at(bond_tool);
    m_cycle_on_click = bond_tool == BondTool::SINGLE;
}

void DrawBondSceneTool::onAtomClicked(const RDKit::Atom* const atom)
{
    AbstractDrawBondSceneTool::onAtomClicked(atom);
    m_last_bond_clicked = nullptr;
}

void DrawBondSceneTool::onBondClicked(const RDKit::Bond* const bond)
{
    if (m_cycle_on_click && m_last_bond_clicked == bond) {
        cycleBond(bond);
    } else {
        AbstractDrawBondSceneTool::onBondClicked(bond);
    }
    m_last_bond_clicked = bond;
}

void DrawBondSceneTool::onEmptySpaceClicked(const RDGeom::Point3D& pos)
{
    AbstractDrawBondSceneTool::onEmptySpaceClicked(pos);
    m_last_bond_clicked = nullptr;
}

bool DrawBondSceneTool::bondMatches(const RDKit::Bond* const bond)
{
    return (!bond->hasQuery() && bond->getBondType() == m_bond_type &&
            bond->getBondDir() == m_bond_dir);
}

void DrawBondSceneTool::mutateBond(const RDKit::Bond* const bond)
{
    m_mol_model->mutateBond(bond, m_bond_type, m_bond_dir);
}

void DrawBondSceneTool::addAtom(const RDGeom::Point3D& pos,
                                const RDKit::Atom* const bound_to)
{
    m_mol_model->addAtom(Element::C, pos, bound_to, m_bond_type, m_bond_dir);
}

void DrawBondSceneTool::addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                                         const RDGeom::Point3D& pos2)
{
    m_mol_model->addAtomChain(Element::C, {pos1, pos2}, nullptr, m_bond_type,
                              m_bond_dir);
}

void DrawBondSceneTool::addBond(const RDKit::Atom* const start_atom,
                                const RDKit::Atom* const end_atom)
{
    m_mol_model->addBond(start_atom, end_atom, m_bond_type, m_bond_dir);
}

DrawBondQuerySceneTool::DrawBondQuerySceneTool(BondTool bond_tool, Scene* scene,
                                               MolModel* mol_model) :
    AbstractDrawBondSceneTool(scene, mol_model)
{
    std::tie(m_query_type, m_query_func) = BOND_TOOL_QUERY_MAP.at(bond_tool);
}

std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>
DrawBondQuerySceneTool::getQuery()
{
    auto query =
        std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>(m_query_func());
    query->setTypeLabel(m_query_type);
    return query;
}

bool DrawBondQuerySceneTool::bondMatches(const RDKit::Bond* const bond)
{
    return bond->hasQuery() && bond->getQuery()->getTypeLabel() == m_query_type;
}

void DrawBondQuerySceneTool::mutateBond(const RDKit::Bond* const bond)
{
    m_mol_model->mutateBond(bond, getQuery());
}

void DrawBondQuerySceneTool::addAtom(const RDGeom::Point3D& pos,
                                     const RDKit::Atom* const bound_to)
{
    m_mol_model->addAtom(Element::C, pos, bound_to, getQuery());
}

void DrawBondQuerySceneTool::addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                                              const RDGeom::Point3D& pos2)
{
    m_mol_model->addAtomChain(Element::C, {pos1, pos2}, nullptr, getQuery());
}

void DrawBondQuerySceneTool::addBond(const RDKit::Atom* const start_atom,
                                     const RDKit::Atom* const end_atom)
{
    m_mol_model->addBond(start_atom, end_atom, getQuery());
}

} // namespace sketcher
} // namespace schrodinger
