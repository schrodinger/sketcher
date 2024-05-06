#pragma once
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

namespace RDKit
{
class Atom;
class QueryAtom;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

enum class QueryType {
    ALLOWED_LIST,
    NOT_ALLOWED_LIST,
    WILDCARD,
    SPECIFIC_ELEMENT,
    RGROUP,
};

enum class EnhancedStereoType { ABS, AND, OR };

enum class QueryAromaticity {
    ANY,
    AROMATIC,
    ALIPHATIC,
};

enum class QueryCount {
    ANY,
    POSITIVE,
    EXACTLY,
};

/**
 * A data class representing atom properties that can be specified in the Edit
 * Atom Properties dialog
 */
struct SKETCHER_API AbstractAtomProperties {
    Element element = Element::C;
    std::optional<unsigned int> isotope = std::nullopt;
    std::optional<EnhancedStereoType> enhanced_stereo_type = std::nullopt;
    unsigned int enhanced_stereo_group_id = 0;

    virtual ~AbstractAtomProperties() = default;
    virtual bool isQuery() const = 0;
    bool operator==(const AbstractAtomProperties& other) const;
};

/**
 * Atom properties that can be specified in the atom (i.e. non-query) page of
 * the Edit Atom Properties dialog
 */
struct SKETCHER_API AtomProperties : AbstractAtomProperties {
    int charge = 0;
    unsigned int unpaired_electrons = 0;

    bool isQuery() const
    {
        return false;
    };
};

/**
 * Atom query properties that can be specified in the query page of the Edit
 * Atom Properties dialog
 */
struct SKETCHER_API AtomQueryProperties : AbstractAtomProperties {
    std::optional<int> charge = std::nullopt;
    std::optional<unsigned int> unpaired_electrons = std::nullopt;
    QueryType query_type = QueryType::ALLOWED_LIST;
    std::vector<Element> allowed_list;
    AtomQuery wildcard = AtomQuery::A;
    unsigned int r_group = 1;
    std::optional<unsigned int> total_h = std::nullopt;
    std::optional<unsigned int> num_connections = std::nullopt;
    QueryAromaticity aromaticity = QueryAromaticity::ANY;
    QueryCount ring_count_type = QueryCount::ANY;
    unsigned int ring_count_exact_val = 0;
    QueryCount ring_bond_count_type = QueryCount::ANY;
    unsigned int ring_bond_count_exact_val = 0;
    std::optional<unsigned int> smallest_ring_size = std::nullopt;
    std::string smarts_query;

    bool isQuery() const
    {
        return true;
    };
};

/**
 * @return the properties for the specified atom
 */
SKETCHER_API std::shared_ptr<AbstractAtomProperties>
read_properties_from_atom(const RDKit::Atom* const atom);

/**
 * @return a new atom with the specified properties
 */
SKETCHER_API RDKit::Atom create_atom_with_properties(
    const std::shared_ptr<AbstractAtomProperties> properties);

} // namespace sketcher
} // namespace schrodinger
