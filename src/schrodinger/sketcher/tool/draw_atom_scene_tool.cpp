#include "schrodinger/sketcher/tool/draw_atom_scene_tool.h"

#include <unordered_map>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/QueryOps.h>

#include <QPainter>

#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"

namespace schrodinger
{
namespace sketcher
{

AbstractDrawAtomSceneTool::AbstractDrawAtomSceneTool(const Fonts& fonts,
                                                     Scene* scene,
                                                     MolModel* mol_model) :
    AbstractDrawSceneTool(scene, mol_model),
    m_fonts(&fonts)
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

DrawElementSceneTool::DrawElementSceneTool(const Element element,
                                           const Fonts& fonts, Scene* scene,
                                           MolModel* mol_model) :
    AbstractDrawAtomSceneTool(fonts, scene, mol_model),
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
    m_mol_model->mutateAtoms({atom}, m_element);
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

QPixmap DrawElementSceneTool::createDefaultCursorPixmap() const
{
    auto elem_name =
        atomic_number_to_symbol(static_cast<unsigned int>(m_element));
    auto qelem_name = QString::fromStdString(elem_name);
    return render_text_to_pixmap(qelem_name, m_fonts->m_cursor_hint_font,
                                 CURSOR_HINT_COLOR);
}

DrawAtomQuerySceneTool::DrawAtomQuerySceneTool(AtomQuery atom_query,
                                               const Fonts& fonts, Scene* scene,
                                               MolModel* mol_model) :
    AbstractDrawAtomSceneTool(fonts, scene, mol_model)
{
    m_atom_query = atom_query;
    m_query_type = getQuery()->getTypeLabel();
}

std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY>
DrawAtomQuerySceneTool::getQuery()
{
    auto query_func = ATOM_TOOL_QUERY_MAP.at(m_atom_query);
    return std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY>(query_func());
}

bool DrawAtomQuerySceneTool::shouldDrawBondForClickOnAtom(
    const RDKit::Atom* const atom) const
{
    return (atom->hasQuery() && atom->getQueryType() == m_query_type);
}

void DrawAtomQuerySceneTool::mutateAtom(const RDKit::Atom* const atom)
{
    m_mol_model->mutateAtoms({atom}, m_atom_query);
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

QPixmap DrawAtomQuerySceneTool::createDefaultCursorPixmap() const
{
    QFont italics = m_fonts->m_cursor_hint_font;
    italics.setItalic(true);
    return render_text_to_pixmap(QString::fromStdString(m_query_type), italics,
                                 CURSOR_HINT_COLOR);
}

} // namespace sketcher
} // namespace schrodinger
