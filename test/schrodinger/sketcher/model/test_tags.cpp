
#define BOOST_TEST_MODULE tags

#include <unordered_map>

#include "../test_common.h"
#include "schrodinger/sketcher/model/tags.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Test basic functionality for the tag classes
 */
BOOST_AUTO_TEST_CASE(test_tags)
{
    auto atom_tag = AtomTag(0);
    BOOST_TEST(atom_tag == AtomTag(0));
    BOOST_TEST(atom_tag != AtomTag(1));
    BOOST_TEST(atom_tag == 0);
    BOOST_TEST(atom_tag != 1);
    BOOST_TEST(atom_tag++ == AtomTag(0));
    BOOST_TEST(atom_tag == AtomTag(1));
    BOOST_TEST(++atom_tag == AtomTag(2));
    BOOST_TEST(atom_tag == AtomTag(2));
    atom_tag += 3;
    BOOST_TEST(atom_tag == AtomTag(5));
    BOOST_TEST(atom_tag.toString() == "AtomTag(5)");
    BOOST_TEST(static_cast<int>(atom_tag) == 5);
}

/**
 * Make sure that tags work properly with maps (which implies that they will
 * also work properly with sets)
 */
BOOST_AUTO_TEST_CASE(test_tag_hashing)
{
    std::unordered_map<AtomTag, int> map;
    map[AtomTag(2)] = 12;
    map[AtomTag(5)] = 15;
    BOOST_TEST(map[AtomTag(2)] == 12);
    BOOST_TEST(map[AtomTag(5)] == 15);
    BOOST_TEST(map.count(AtomTag(5)) == 1);
    BOOST_TEST(map.count(AtomTag(6)) == 0);
}

} // namespace sketcher
} // namespace schrodinger