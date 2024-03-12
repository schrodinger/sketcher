#pragma once

#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/tool/draw_atom_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A base class for scene tools that create R-group atoms
 */
class SKETCHER_API AbstractDrawRGroupSceneTool
    : public AbstractDrawAtomSceneTool
{
  public:
    AbstractDrawRGroupSceneTool(const Fonts& fonts, Scene* scene,
                                MolModel* mol_model);

  protected:
    /**
     * Get R-group numbers to use for new R-group atoms
     * @param how_many How many R-group numbers to return
     */
    virtual std::vector<unsigned int>
    getRGroupNums(const size_t how_many) const = 0;

    // reimplemented AbstractDrawSceneTool methods
    void mutateAtom(const RDKit::Atom* const atom) override;
    void addAtom(const RDGeom::Point3D& pos,
                 const RDKit::Atom* const bound_to = nullptr) override;
    void addTwoBoundAtoms(const RDGeom::Point3D& pos1,
                          const RDGeom::Point3D& pos2) override;
};

/**
 * A scene tool for drawing R-group atoms with a specific R-group number
 */
class SKETCHER_API DrawRGroupSceneTool : public AbstractDrawRGroupSceneTool
{
  public:
    DrawRGroupSceneTool(const int r_group_num, const Fonts& fonts, Scene* scene,
                        MolModel* mol_model);

  protected:
    unsigned int m_r_group_num;

    // reimplemented AbstractDrawRGroupSceneTool methods
    std::vector<unsigned int>
    getRGroupNums(const size_t how_many) const override;

    // reimplemented AbstractDrawSceneTool methods
    bool
    shouldDrawBondForClickOnAtom(const RDKit::Atom* const atom) const override;

    // reimplemented AbstractSceneTool methods
    QPixmap createDefaultCursorPixmap() const override;
};

/**
 * A scene tool for drawing R-group atoms that uses the smallest unused R-group
 * number
 */
class SKETCHER_API DrawIncrementingRGroupSceneTool
    : public AbstractDrawRGroupSceneTool
{
  public:
    DrawIncrementingRGroupSceneTool(const Fonts& fonts, Scene* scene,
                                    MolModel* mol_model);

  protected:
    // reimplemented AbstractDrawRGroupSceneTool methods
    std::vector<unsigned int>
    getRGroupNums(const size_t how_many) const override;

    // reimplemented AbstractDrawSceneTool methods
    bool
    shouldDrawBondForClickOnAtom(const RDKit::Atom* const atom) const override;
    QPixmap createDefaultCursorPixmap() const override;
};

} // namespace sketcher
} // namespace schrodinger
