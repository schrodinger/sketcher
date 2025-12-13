#define BOOST_TEST_MODULE Test_Sketcher

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/rdkit/coord_utils.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/tool/move_rotate_scene_tool.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

BOOST_AUTO_TEST_CASE(test_get_overlapping_atom_idxs)
{
    auto mol = rdkit_extensions::to_rdkit("C.C.C");
    sketcher::update_2d_coordinates(*mol);
    auto& conf = mol->getConformer();
    conf.setAtomPos(0, {-MAX_DIST_FOR_DRAG_MERGE / 2, 0, 0});
    conf.setAtomPos(1, {MAX_DIST_FOR_DRAG_MERGE / 2, 0, 0});
    conf.setAtomPos(2, {-MAX_DIST_FOR_DRAG_MERGE / 4, 0, 0});
    std::unordered_set<const RDKit::Atom*> atoms_to_move = {
        mol->getAtomWithIdx(2)};

    // atom 2 is overlapping with both atoms 0 and 1, but it's closer to atom 0
    auto overlapping = get_overlapping_atom_idxs(&*mol, atoms_to_move);
    std::vector<std::pair<unsigned int, unsigned int>> exp_overlapping = {
        {2, 0}};
    BOOST_TEST(overlapping == exp_overlapping);

    // atom 2 is now closer to atom 1
    conf.setAtomPos(2, {MAX_DIST_FOR_DRAG_MERGE / 4, 0, 0});
    overlapping = get_overlapping_atom_idxs(&*mol, atoms_to_move);
    exp_overlapping = {{2, 1}};
    BOOST_TEST(overlapping == exp_overlapping);

    // now atom 2 isn't overlapping with any atoms
    conf.setAtomPos(2, {2 * MAX_DIST_FOR_DRAG_MERGE, 0, 0});
    overlapping = get_overlapping_atom_idxs(&*mol, atoms_to_move);
    exp_overlapping = {};
    BOOST_TEST(overlapping == exp_overlapping);
}

} // namespace sketcher
} // namespace schrodinger
