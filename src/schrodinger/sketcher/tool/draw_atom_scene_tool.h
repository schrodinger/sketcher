#pragma once

#include <functional>
#include <memory>
#include <string>

#include <rdkit/GraphMol/QueryAtom.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/tool/abstract_draw_atom_bond_scene_tool.h"

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace RDGeom
{
class Point3D;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * A base class for scene tools that draw atoms or query atoms
 */
class SKETCHER_API AbstractDrawAtomSceneTool : public AbstractDrawSceneTool
{
  public:
    AbstractDrawAtomSceneTool(const Fonts& fonts, Scene* scene,
                              MolModel* mol_model);

  protected:
    const Fonts* m_fonts;

    // reimplemented AbstractDrawSceneTool methods
    void onBondClicked(const RDKit::Bond* const bond) override;
    void onEmptySpaceClicked(const RDGeom::Point3D& pos) override;
    void addBond(const RDKit::Atom* const start_atom,
                 const RDKit::Atom* const end_atom) override;
};

/**
 * A scene tool for drawing elements (i.e. non-query atoms)
 */
class SKETCHER_API DrawElementSceneTool : public AbstractDrawAtomSceneTool
{
  public:
    DrawElementSceneTool(const Element element, const Fonts& fonts,
                         Scene* scene, MolModel* mol_model);

  protected:
    Element m_element;

    // reimplemented AbstractDrawSceneTool methods
    bool
    shouldDrawBondForClickOnAtom(const RDKit::Atom* const atom) const override;
    void mutateAtom(const RDKit::Atom* const atom) override;
    void addAtom(const RDGeom::Point3D& pos,
                 const RDKit::Atom* const bound_to = nullptr) override;
    void addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                          const RDGeom::Point3D& pos2) override;

    // reimplemented AbstractSceneTool method
    QPixmap getCursorPixmap() const override;
};

/**
 * A scene tool for drawing query atoms
 */
class SKETCHER_API DrawAtomQuerySceneTool : public AbstractDrawAtomSceneTool
{
  public:
    DrawAtomQuerySceneTool(AtomQuery atom_query, const Fonts& fonts,
                           Scene* scene, MolModel* mol_model);

  protected:
    /**
     * A function that returns an RDKit Query object containing the appropriate
     * atom query
     */
    std::function<RDKit::QueryAtom::QUERYATOM_QUERY*()> m_query_func;
    std::string m_query_type;

    /**
     * @return the query to use for atoms
     */
    std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> getQuery();

    // reimplemented AbstractDrawSceneTool methods
    bool
    shouldDrawBondForClickOnAtom(const RDKit::Atom* const atom) const override;
    void mutateAtom(const RDKit::Atom* const atom) override;
    void addAtom(const RDGeom::Point3D& pos,
                 const RDKit::Atom* const bound_to = nullptr) override;
    void addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                          const RDGeom::Point3D& pos2) override;
    QPixmap getCursorPixmap() const override;
};

} // namespace sketcher
} // namespace schrodinger
