#include "schrodinger/sketcher/model/tags.h"

#include <functional>
#include <unordered_map>
#include <ostream>

namespace schrodinger
{
namespace sketcher
{

template <TagType tag_type> Tag<tag_type>::Tag(const int val) noexcept :
    m_val(val)
{
}

template <TagType tag_type> Tag<tag_type>::operator const int&() const
{
    return m_val;
}

template <TagType tag_type> Tag<tag_type>::operator int&()
{
    return m_val;
}

template <TagType tag_type>
bool Tag<tag_type>::operator==(const Tag& other) const
{
    return m_val == other.m_val;
}

template <TagType tag_type>
bool Tag<tag_type>::operator<(const Tag& other) const
{
    return m_val < other.m_val;
}

template <TagType tag_type> Tag<tag_type>& Tag<tag_type>::operator++()
{
    ++m_val;
    return *this;
}

template <TagType tag_type>
Tag<tag_type> Tag<tag_type>::operator++(int increment_by)
{
    if (increment_by == 0) {
        // tag++ will call this method with a 0
        increment_by = 1;
    }
    Tag<tag_type> before = *this;
    m_val += increment_by;
    return before;
}

template <> std::string Tag<TagType::ATOM>::toString() const
{
    return "AtomTag(" + std::to_string(m_val) + ")";
}

template <> std::string Tag<TagType::BOND>::toString() const
{
    return "BondTag(" + std::to_string(m_val) + ")";
}

template <> std::string Tag<TagType::S_GROUP>::toString() const
{
    return "SGroupTag(" + std::to_string(m_val) + ")";
}

template <> std::string Tag<TagType::NON_MOLECULAR>::toString() const
{
    return "NonMolecularTag(" + std::to_string(m_val) + ")";
}

template <TagType tag_type>
std::ostream& operator<<(std::ostream& out, const Tag<tag_type>& tag)
{
    return out << tag.toString();
}

// explicitly instantiate all templates since we define them in the
// implementation file rather than the header

template class Tag<TagType::ATOM>;
template class Tag<TagType::BOND>;
template class Tag<TagType::S_GROUP>;
template class Tag<TagType::NON_MOLECULAR>;

template SKETCHER_API std::ostream& operator<<(std::ostream& out,
                                               const AtomTag& tag);
template SKETCHER_API std::ostream& operator<<(std::ostream& out,
                                               const BondTag& tag);
template SKETCHER_API std::ostream& operator<<(std::ostream& out,
                                               const SGroupTag& tag);
template SKETCHER_API std::ostream& operator<<(std::ostream& out,
                                               const NonMolecularTag& tag);

} // namespace sketcher
} // namespace schrodinger

using schrodinger::sketcher::Tag;
using schrodinger::sketcher::TagType;

template <TagType tag_type> std::size_t
std::hash<Tag<tag_type>>::operator()(const Tag<tag_type>& tag) const noexcept
{
    return std::hash<int>{}(tag);
}

template struct SKETCHER_API std::hash<Tag<TagType::ATOM>>;
template struct SKETCHER_API std::hash<Tag<TagType::BOND>>;
template struct SKETCHER_API std::hash<Tag<TagType::S_GROUP>>;
template struct SKETCHER_API std::hash<Tag<TagType::NON_MOLECULAR>>;
