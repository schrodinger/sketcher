#include "schrodinger/sketcher/tool/draw_bond_scene_tool.h"

#include <tuple>
#include <unordered_map>

#include <rdkit/GraphMol/QueryBond.h>
#include <rdkit/GraphMol/QueryOps.h>

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
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
    m_bond_tool = bond_tool;
    std::tie(m_bond_type, m_bond_dir, m_flippable, m_cursor_hint_path) =
        BOND_TOOL_BOND_MAP.at(bond_tool);
    m_cycle_on_click = bond_tool == BondTool::SINGLE;
}

void DrawBondSceneTool::onAtomClicked(const RDKit::Atom* const atom)
{
    AbstractDrawBondSceneTool::onAtomClicked(atom);
    m_last_clicked_bond_idx = -1;
}

void DrawBondSceneTool::onBondClicked(const RDKit::Bond* const bond)
{
    // save the index of the bond before we call the command and it gets
    // invalidated
    auto bond_idx = bond->getIdx();
    if (m_cycle_on_click &&
        m_last_clicked_bond_idx == static_cast<int>(bond->getIdx())) {
        cycleBond(bond);
    } else {
        AbstractDrawBondSceneTool::onBondClicked(bond);
    }
    m_last_clicked_bond_idx = bond_idx;
}

void DrawBondSceneTool::onEmptySpaceClicked(const RDGeom::Point3D& pos)
{
    AbstractDrawBondSceneTool::onEmptySpaceClicked(pos);
    m_last_clicked_bond_idx = -1;
}

bool DrawBondSceneTool::bondMatches(const RDKit::Bond* const bond)
{
    return (!bond->hasQuery() && bond->getBondType() == m_bond_type &&
            bond->getBondDir() == m_bond_dir);
}

void DrawBondSceneTool::mutateBond(const RDKit::Bond* const bond)
{
    m_mol_model->mutateBonds({bond}, m_bond_tool);
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

QPixmap DrawBondSceneTool::createDefaultCursorPixmap() const
{
    return cursor_hint_from_svg(m_cursor_hint_path);
}

DrawBondQuerySceneTool::DrawBondQuerySceneTool(BondTool bond_tool,
                                               const Fonts& fonts, Scene* scene,
                                               MolModel* mol_model) :
    AbstractDrawBondSceneTool(scene, mol_model),
    m_fonts(&fonts)
{
    m_bond_tool = bond_tool;
    m_query_type = get_label_for_bond_query(getQueryAndBondType().first.get());
}

std::pair<std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>,
          RDKit::Bond::BondType>
DrawBondQuerySceneTool::getQueryAndBondType()
{
    auto [query_func, bond_type] = BOND_TOOL_QUERY_MAP.at(m_bond_tool);
    auto query =
        std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>(query_func());
    return {query, bond_type};
}

bool DrawBondQuerySceneTool::bondMatches(const RDKit::Bond* const bond)
{
    auto [bond_type, query_type] = get_bond_type_and_query_label(bond);
    return query_type == m_query_type;
}

void DrawBondQuerySceneTool::mutateBond(const RDKit::Bond* const bond)
{
    m_mol_model->mutateBonds({bond}, m_bond_tool);
}

void DrawBondQuerySceneTool::addAtom(const RDGeom::Point3D& pos,
                                     const RDKit::Atom* const bound_to)
{
    auto [query, bond_type] = getQueryAndBondType();
    m_mol_model->addAtom(Element::C, pos, bound_to, query, bond_type);
}

void DrawBondQuerySceneTool::addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                                              const RDGeom::Point3D& pos2)
{
    auto [query, bond_type] = getQueryAndBondType();
    m_mol_model->addAtomChain(Element::C, {pos1, pos2}, nullptr, query,
                              bond_type);
}

void DrawBondQuerySceneTool::addBond(const RDKit::Atom* const start_atom,
                                     const RDKit::Atom* const end_atom)
{
    auto [query, bond_type] = getQueryAndBondType();
    m_mol_model->addBond(start_atom, end_atom, query, bond_type);
}

QPixmap DrawBondQuerySceneTool::createDefaultCursorPixmap() const
{
    auto text = QString::fromStdString(m_query_type);
    return render_text_to_pixmap(text, m_fonts->m_cursor_hint_font,
                                 CURSOR_HINT_COLOR);
}

} // namespace sketcher
} // namespace schrodinger
