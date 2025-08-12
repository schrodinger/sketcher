/**
 * Definitions for tag types using in MolModel: AtomTag, BondTag, SGroupTag, and
 * NonMolecularTag.  These tags work like ints, but with additional type safety
 * so that, e.g., AtomTags can't be passed to a function that expects BondTags.
 */

#pragma once

#include <string>

#include <boost/operators.hpp>
#include <boost/serialization/strong_typedef.hpp>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * An enum used to distinguish the different types of tags.  This enum should
 * never need to be used outside of this module.  Instead, use the typedefs
 * defined below (AtomTag, BondTag, SGroupTag, and NonMolecularTag) to access
 * the tag classes.
 */
enum class TagType { ATOM, BOND, S_GROUP, NON_MOLECULAR };

template <TagType tag_type> class SKETCHER_API Tag
    : public boost::totally_ordered<Tag<tag_type>>
{
  public:
    Tag() noexcept = default;
    explicit Tag(const int val) noexcept;

    operator const int&() const;
    operator int&();

    bool operator==(const Tag& other) const;
    bool operator<(const Tag& other) const;
    Tag& operator++();
    Tag operator++(int increment_by);

    /**
     * @return the tag as a string, e.g. "AtomTag(5)".  This functionality is
     * primarily intended for use when debugging.
     */
    std::string toString() const;

  private:
    int m_val = 0;
};

typedef Tag<TagType::ATOM> AtomTag;
typedef Tag<TagType::BOND> BondTag;
typedef Tag<TagType::S_GROUP> SGroupTag;
typedef Tag<TagType::NON_MOLECULAR> NonMolecularTag;

/**
 * An output stream operator for the tag classes.  This functionality is
 * primarily intended for use when debugging.
 */
template <TagType tag_type> SKETCHER_API std::ostream&
operator<<(std::ostream& out, const Tag<tag_type>& tag);

} // namespace sketcher
} // namespace schrodinger

// define hashing so we can make sets and maps of tags
template <schrodinger::sketcher::TagType tag_type>
struct std::hash<schrodinger::sketcher::Tag<tag_type>> {
    std::size_t
    operator()(const schrodinger::sketcher::Tag<tag_type>& tag) const noexcept;
};
