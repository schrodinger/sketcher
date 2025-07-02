#pragma once

#include <boost/dynamic_bitset.hpp>
#include <iterator>

namespace schrodinger
{

/**
 * An iterator which walks through the "on" positions in a
 * boost::dynamic_bitset
 */
class dynamic_bitset_on_bits_iterator
{
  public:
    using iterator_category = std::forward_iterator_tag;
    using difference_type = std::ptrdiff_t;
    using value_type = size_t;
    using pointer = size_t*;
    using reference = size_t&;

    dynamic_bitset_on_bits_iterator() :
        dynamic_bitset_on_bits_iterator(nullptr, boost::dynamic_bitset<>::npos)
    {
    }

    dynamic_bitset_on_bits_iterator(const boost::dynamic_bitset<>* bitset) :
        dynamic_bitset_on_bits_iterator(bitset, bitset->find_first())
    {
    }

    dynamic_bitset_on_bits_iterator(const boost::dynamic_bitset<>* bitset,
                                    size_t position) :
        d_bitset_p(bitset),
        d_position(position)
    {
    }

    /// Copy constructor
    dynamic_bitset_on_bits_iterator(const dynamic_bitset_on_bits_iterator& it) :
        d_bitset_p(it.d_bitset_p),
        d_position(it.d_position)
    {
    }

    dynamic_bitset_on_bits_iterator& operator++()
    {
        d_position = d_bitset_p->find_next(d_position);
        return *this;
    }

    value_type operator*() const
    {
        return d_position;
    }

    bool operator==(const dynamic_bitset_on_bits_iterator& it) const
    {
        return d_position == it.d_position;
    }

    bool operator!=(const dynamic_bitset_on_bits_iterator& it) const
    {
        return !(*this == it);
    }

  private:
    const boost::dynamic_bitset<>* d_bitset_p;
    size_t d_position;
};

/**
 * A helper class to iterate over the "on" positions in a boost::dynamic_bitset.
 * This mostly serve to provide `begin()` and `end()` members for `for loops`
 * and other standard algorithms.
 */
class dynamic_bitset_on_bits_wrapper
{
  public:
    dynamic_bitset_on_bits_wrapper() = delete;
    dynamic_bitset_on_bits_wrapper(const boost::dynamic_bitset<>& bitset) :
        d_bitset_p(&bitset)
    {
    }

    using iterator = dynamic_bitset_on_bits_iterator;

    iterator begin() const
    {
        return iterator(d_bitset_p);
    }

    iterator end() const
    {
        return iterator(d_bitset_p, boost::dynamic_bitset<>::npos);
    }

  private:
    const boost::dynamic_bitset<>* d_bitset_p;
};

} // namespace schrodinger
