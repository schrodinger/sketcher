
#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/rdkit/s_group_constants.h"
#include "schrodinger/sketcher/dialog/bracket_subgroup_dialog.h"
#include "schrodinger/sketcher/model/mol_model.h"
#include "schrodinger/sketcher/ui/ui_bracket_subgroup_dialog.h"

#include <QUndoStack>

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::SubgroupType)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::RepeatPattern)

namespace schrodinger
{
namespace sketcher
{

class TestBracketSubgroupDialog : public BracketSubgroupDialog
{
  public:
    TestBracketSubgroupDialog() : BracketSubgroupDialog(&m_test_mol_model){};
    QUndoStack m_undo_stack;
    MolModel m_test_mol_model = MolModel(&m_undo_stack);

    using BracketSubgroupDialog::getPolymerLabel;
    using BracketSubgroupDialog::getRepeatPattern;
    using BracketSubgroupDialog::getSubgroupType;
};

/**
 * Verify that the substance group type can be assigned.
 */
BOOST_AUTO_TEST_CASE(test_setSubgroupType)
{
    TestBracketSubgroupDialog dialog;

    // Verify the default behavior
    BOOST_TEST(static_cast<int>(dialog.getSubgroupType()) ==
               static_cast<int>(SubgroupType::SRU_POLYMER));

    dialog.setSubgroupType(SubgroupType::COPOLYMER);
    BOOST_TEST(static_cast<int>(dialog.getSubgroupType()) ==
               static_cast<int>(SubgroupType::COPOLYMER));

    // If the assigned type is unrecognized, use the default
    dialog.setSubgroupType(SubgroupType::OTHER);
    BOOST_TEST(static_cast<int>(dialog.getSubgroupType()) ==
               static_cast<int>(SubgroupType::SRU_POLYMER));
}

/**
 * Verify that the repeating pattern type can be assigned.
 */
BOOST_AUTO_TEST_CASE(test_setRepeatPattern)
{
    TestBracketSubgroupDialog dialog;

    // Verify the default behavior
    BOOST_TEST(static_cast<int>(dialog.getRepeatPattern()) ==
               static_cast<int>(RepeatPattern::HEAD_TO_TAIL));

    dialog.setRepeatPattern(RepeatPattern::HEAD_TO_TAIL);
    BOOST_TEST(static_cast<int>(dialog.getRepeatPattern()) ==
               static_cast<int>(RepeatPattern::HEAD_TO_TAIL));
    dialog.setRepeatPattern(RepeatPattern::HEAD_TO_HEAD);
    BOOST_TEST(static_cast<int>(dialog.getRepeatPattern()) ==
               static_cast<int>(RepeatPattern::HEAD_TO_HEAD));
    dialog.setRepeatPattern(RepeatPattern::EITHER_UNKNOWN);
    BOOST_TEST(static_cast<int>(dialog.getRepeatPattern()) ==
               static_cast<int>(RepeatPattern::EITHER_UNKNOWN));
}

/**
 * Verify that the polymer label can be assigned.
 */
BOOST_AUTO_TEST_CASE(test_setPolymerLabel)
{
    TestBracketSubgroupDialog dialog;

    // Verify the default behavior
    BOOST_TEST(dialog.getPolymerLabel().toStdString() == "");

    // If assigning a copolymer, the polymer label text gets a default value
    dialog.setSubgroupType(SubgroupType::COPOLYMER);
    BOOST_TEST(dialog.getPolymerLabel().toStdString() == "co");

    // This value cannot be changed until the subgroup type changes back to one
    // that allows polymer label editing
    dialog.setPolymerLabel("n");
    BOOST_TEST(dialog.getPolymerLabel().toStdString() == "co");

    dialog.setSubgroupType(SubgroupType::SRU_POLYMER);
    BOOST_TEST(dialog.getPolymerLabel().toStdString() == "");

    dialog.setPolymerLabel("n");
    BOOST_TEST(dialog.getPolymerLabel().toStdString() == "n");

    // We also use the validator to prevent certain values from being assigned
    dialog.setPolymerLabel("abc");
    BOOST_TEST(dialog.getPolymerLabel().toStdString() == "n");
}

/**
 * Make sure that setSubgroup loads the subgroup's current settings into the
 * dialog
 */
BOOST_AUTO_TEST_CASE(test_setSubgroup)
{
    TestBracketSubgroupDialog dialog;
    auto mol = rdkit_extensions::to_rdkit("CCCCCCCCC |Sg:n:3,4,5,6:7:hh:::|");
    auto& s_group = RDKit::getSubstanceGroups(*mol)[0];
    dialog.setSubgroup(&s_group);
    BOOST_TEST(dialog.getSubgroupType() == SubgroupType::SRU_POLYMER);
    BOOST_TEST(dialog.getPolymerLabel().toStdString() == "7");
    BOOST_TEST(dialog.getRepeatPattern() == RepeatPattern::HEAD_TO_HEAD);

    auto mol2 = rdkit_extensions::to_rdkit("CCCCCCCCC |Sg:alt:1,2,3:co:eu:::|");
    auto& s_group2 = RDKit::getSubstanceGroups(*mol2)[0];
    dialog.setSubgroup(&s_group2);
    BOOST_TEST(dialog.getSubgroupType() == SubgroupType::COPOLYMER);
    BOOST_TEST(dialog.getPolymerLabel().toStdString() == "co");
    BOOST_TEST(dialog.getRepeatPattern() == RepeatPattern::EITHER_UNKNOWN);
}

} // namespace sketcher
} // namespace schrodinger
