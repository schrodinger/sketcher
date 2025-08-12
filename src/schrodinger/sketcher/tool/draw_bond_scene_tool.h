#pragma once

#include <functional>
#include <memory>
#include <string>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/QueryBond.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/tool/abstract_draw_atom_bond_scene_tool.h"

namespace RDGeom
{
class Point3D;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * A base class for scene tools that draw bonds or query bonds
 */
class SKETCHER_API AbstractDrawBondSceneTool : public AbstractDrawSceneTool
{
  public:
    AbstractDrawBondSceneTool(Scene* scene, MolModel* mol_model);

  protected:
    /**
     * Does the bond change meaningfully if we swap the start and end atoms?
     */
    bool m_flippable = false;

    /**
     * Should the scene tool cycle bonds (i.e. single -> double -> triple ->
     * single...) when clicking on the same bond multiple times in a row.
     * Subclasses that use this variable (e.g. DrawBondSceneTool) are
     * responsible for updating it.
     */
    bool m_cycle_on_click;

    // overriden AbstractDrawSceneTool methods
    bool
    shouldDrawBondForClickOnAtom(const RDKit::Atom* const atom) const override;
    virtual void onBondClicked(const RDKit::Bond* const bond) override;
    virtual void onEmptySpaceClicked(const RDGeom::Point3D& pos) override;

    /**
     * @return Whether the specified bond is identical (same BondType and same
     * BondDir, or identical query) to a bond that this tool would create.
     */
    virtual bool bondMatches(const RDKit::Bond* const bond) = 0;

    /**
     * Mutate the specified bond to match the bonds that this tool would create.
     */
    virtual void mutateBond(const RDKit::Bond* const bond) = 0;
};

/**
 * A scene tool for drawing bonds other than query bonds
 */
class SKETCHER_API DrawBondSceneTool : public AbstractDrawBondSceneTool
{
  public:
    DrawBondSceneTool(BondTool bond_tool, Scene* scene, MolModel* mol_model);

  protected:
    BondTool m_bond_tool;
    RDKit::Bond::BondType m_bond_type;
    RDKit::Bond::BondDir m_bond_dir;
    int m_last_clicked_bond_idx = -1;
    QString m_cursor_hint_path;

    // overriden AbstractDrawBondSceneTool methods
    bool bondMatches(const RDKit::Bond* const bond) override;
    void mutateBond(const RDKit::Bond* const bond) override;

    // overriden AbstractDrawSceneTool methods
    void onAtomClicked(const RDKit::Atom* const atom) override;
    void onBondClicked(const RDKit::Bond* const bond) override;
    void onEmptySpaceClicked(const RDGeom::Point3D& pos) override;
    void addAtom(const RDGeom::Point3D& pos,
                 const RDKit::Atom* const bound_to = nullptr) override;
    void addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                          const RDGeom::Point3D& pos2) override;
    void addBond(const RDKit::Atom* const start_atom,
                 const RDKit::Atom* const end_atom) override;

    // overriden AbstractSceneTool method
    QPixmap createDefaultCursorPixmap() const override;
};

/**
 * A scene tool for drawing query bonds
 */
class SKETCHER_API DrawBondQuerySceneTool : public AbstractDrawBondSceneTool
{
  public:
    DrawBondQuerySceneTool(BondTool bond_tool, const Fonts& fonts, Scene* scene,
                           MolModel* mol_model);

  protected:
    /**
     * Sketcher bond type enum
     */
    BondTool m_bond_tool;

    /**
     * The text to use for labeling the drawn queries
     */
    std::string m_query_type;

    const Fonts* m_fonts;

    /**
     * @return the query and bond type to use for bonds
     */
    std::pair<std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>,
              RDKit::Bond::BondType>
    getQueryAndBondType();

    // overriden AbstractDrawBondSceneTool methods
    bool bondMatches(const RDKit::Bond* const bond) override;
    void mutateBond(const RDKit::Bond* const bond) override;

    // overriden AbstractDrawSceneTool methods
    void addAtom(const RDGeom::Point3D& pos,
                 const RDKit::Atom* const bound_to = nullptr) override;
    void addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                          const RDGeom::Point3D& pos2) override;
    void addBond(const RDKit::Atom* const start_atom,
                 const RDKit::Atom* const end_atom) override;

    // overriden AbstractSceneTool methods
    QPixmap createDefaultCursorPixmap() const override;
};

} // namespace sketcher
} // namespace schrodinger
