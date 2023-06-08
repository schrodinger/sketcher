#include "schrodinger/sketcher/tool/draw_atom_scene_tool.h"

#include <unordered_map>

#include <Geometry/point.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryOps.h>

#include "schrodinger/sketcher/model/mol_model.h"

namespace schrodinger
{
namespace sketcher
{

AbstractDrawAtomSceneTool::AbstractDrawAtomSceneTool(Scene* scene,
                                                     MolModel* mol_model) :
    AbstractDrawSceneTool(scene, mol_model)
{
}

void AbstractDrawAtomSceneTool::onBondClicked(const RDKit::Bond* const bond)
{
    cycleBond(bond);
}

void AbstractDrawAtomSceneTool::onEmptySpaceClicked(const RDGeom::Point3D& pos)
{
    addAtom(pos);
}

void AbstractDrawAtomSceneTool::addBond(const RDKit::Atom* const start_atom,
                                        const RDKit::Atom* const end_atom)
{
    m_mol_model->addBond(start_atom, end_atom, RDKit::Bond::BondType::SINGLE);
}

DrawElementSceneTool::DrawElementSceneTool(Element element, Scene* scene,
                                           MolModel* mol_model) :
    AbstractDrawAtomSceneTool(scene, mol_model),
    m_element(element)
{
}

bool DrawElementSceneTool::shouldDrawBondForClickOnAtom(
    const RDKit::Atom* const atom) const
{
    return (!atom->hasQuery() &&
            atom->getAtomicNum() == static_cast<int>(m_element));
}

void DrawElementSceneTool::mutateAtom(const RDKit::Atom* const atom)
{
    m_mol_model->mutateAtom(atom, m_element);
}

void DrawElementSceneTool::addAtom(const RDGeom::Point3D& pos,
                                   const RDKit::Atom* const bound_to)
{
    m_mol_model->addAtom(m_element, pos, bound_to,
                         RDKit::Bond::BondType::SINGLE,
                         RDKit::Bond::BondDir::NONE);
}
void DrawElementSceneTool::addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                                            const RDGeom::Point3D& pos2)
{
    m_mol_model->addAtomChain(m_element, {pos1, pos2}, nullptr,
                              RDKit::Bond::BondType::SINGLE,
                              RDKit::Bond::BondDir::NONE);
}

DrawAtomQuerySceneTool::DrawAtomQuerySceneTool(AtomQuery atom_query,
                                               Scene* scene,
                                               MolModel* mol_model) :
    AbstractDrawAtomSceneTool(scene, mol_model)
{
    const std::unordered_map<
        AtomQuery, std::function<RDKit::QueryAtom::QUERYATOM_QUERY*()>>
        query_type_map = {
            {AtomQuery::A, RDKit::makeAAtomQuery},
            {AtomQuery::AH, RDKit::makeAHAtomQuery},
            {AtomQuery::Q, RDKit::makeQAtomQuery},
            {AtomQuery::QH, RDKit::makeQHAtomQuery},
            {AtomQuery::M, RDKit::makeMAtomQuery},
            {AtomQuery::MH, RDKit::makeMHAtomQuery},
            {AtomQuery::X, RDKit::makeXAtomQuery},
            {AtomQuery::XH, RDKit::makeXHAtomQuery},
        };
    auto query_pair = query_type_map.find(atom_query);
    if (query_pair != query_type_map.end()) {
        m_query_func = query_pair->second;
        auto query = getQuery();
        m_query_type = query->getTypeLabel();
    } else {
        throw std::runtime_error(
            "Invalid AtomQuery passed to DrawAtomQuerySceneTool");
    }
}

std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY>
DrawAtomQuerySceneTool::getQuery()
{
    return std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY>(m_query_func());
}

bool DrawAtomQuerySceneTool::shouldDrawBondForClickOnAtom(
    const RDKit::Atom* const atom) const
{
    return (atom->hasQuery() && atom->getQueryType() == m_query_type);
}

void DrawAtomQuerySceneTool::mutateAtom(const RDKit::Atom* const atom)
{
    m_mol_model->mutateAtom(atom, getQuery());
}

void DrawAtomQuerySceneTool::addAtom(const RDGeom::Point3D& pos,
                                     const RDKit::Atom* const bound_to)
{
    m_mol_model->addAtom(getQuery(), pos, bound_to,
                         RDKit::Bond::BondType::SINGLE,
                         RDKit::Bond::BondDir::NONE);
}

void DrawAtomQuerySceneTool::addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                                              const RDGeom::Point3D& pos2)
{
    m_mol_model->addAtomChain(getQuery(), {pos1, pos2}, nullptr,
                              RDKit::Bond::BondType::SINGLE,
                              RDKit::Bond::BondDir::NONE);
}

} // namespace sketcher
} // namespace schrodinger
