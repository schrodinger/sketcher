#pragma once
#include <iosfwd>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include <rdkit/GraphMol/StereoGroup.h>

#include "schrodinger/rdkit_extensions/stereochemistry.h"
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
    SMARTS,
};

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

// rdkit_extensions::EnhancedStereo is just a std::pair that holds a group type
// enum and id number.  We wrap that class here to add named getters and setters
// in order to make the dialog code more readable
class SKETCHER_API EnhancedStereo : public rdkit_extensions::EnhancedStereo
{
  public:
    using rdkit_extensions::EnhancedStereo::EnhancedStereo;
    EnhancedStereo(rdkit_extensions::EnhancedStereo enh_stereo);
    RDKit::StereoGroupType type() const;
    void setType(RDKit::StereoGroupType);
    unsigned int groupId() const;
    void setGroupId(unsigned int);
};

/**
 * A data class representing atom properties that can be specified in the Edit
 * Atom Properties dialog
 */
struct SKETCHER_API AbstractAtomProperties {
    Element element = Element::C;
    std::optional<unsigned int> isotope = std::nullopt;
    std::optional<EnhancedStereo> enhanced_stereo = std::nullopt;

    virtual ~AbstractAtomProperties() = default;
    virtual bool isQuery() const = 0;
    bool operator==(const AbstractAtomProperties& other) const;
    bool operator!=(const AbstractAtomProperties& other) const;
    friend SKETCHER_API std::ostream&
    operator<<(std::ostream& os, const AbstractAtomProperties& props);
};

/**
 * Print out all contents of an atom properties object. This output is
 * primarily intended for debugging use, not for displaying to the user.
 */
SKETCHER_API std::ostream& operator<<(std::ostream& os,
                                      const AbstractAtomProperties& props);

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
    QueryType query_type = QueryType::SPECIFIC_ELEMENT;
    std::vector<Element> allowed_list;
    AtomQuery wildcard = AtomQuery::A;
    unsigned int r_group = 1;
    std::string smarts_query;
    QueryCount total_h_type = QueryCount::ANY;
    unsigned int total_h_exact_val = 0;
    std::optional<unsigned int> num_connections = std::nullopt;
    QueryAromaticity aromaticity = QueryAromaticity::ANY;
    QueryCount ring_count_type = QueryCount::ANY;
    unsigned int ring_count_exact_val = 0;
    QueryCount ring_bond_count_type = QueryCount::ANY;
    unsigned int ring_bond_count_exact_val = 0;
    std::optional<unsigned int> smallest_ring_size = std::nullopt;

    bool isQuery() const
    {
        return true;
    };

    /**
     * Does this query have any properties other than those that directly relate
     * to the query type (e.g. allowed_list)? Note: stereo properties are
     * excluded from this check since they don't prevent displaying allowed
     * lists with element symbols (SKETCH-2487).
     */
    bool hasPropertiesBeyondQueryType() const;

    /**
     * Does this query have any properties set on the Advanced tab of the Edit
     * Atom Properties dialog?
     */
    bool hasAdvancedProperties() const;
};
/**
 * return a new AtomQueryProperties object with only the advanced properties and
 * default values for the rest
 */
SKETCHER_API std::shared_ptr<AtomQueryProperties>
get_only_advanced_properties(std::shared_ptr<AtomQueryProperties> properties);

/**
 * @return the properties for the specified atom
 */
SKETCHER_API std::shared_ptr<AbstractAtomProperties>
read_properties_from_atom(const RDKit::Atom* const atom);

/**
 * @return the properties for the given query
 */
SKETCHER_API std::shared_ptr<AtomQueryProperties>
read_query(const RDKit::Atom::QUERYATOM_QUERY* const query);

/**
 * @return a new atom with the specified properties along with the settings for
 * the enhanced stereo group that the atom should be added to.  If the atom
 * should not be added to an enhanced stereo group, then that value will be
 * std::nullopt.
 */
SKETCHER_API
std::pair<std::shared_ptr<RDKit::Atom>, std::optional<EnhancedStereo>>
create_atom_with_properties(
    const std::shared_ptr<AbstractAtomProperties> properties);

/**
 * @return the label for the given AtomQuery enum
 */
SKETCHER_API QString get_label_from_atom_query(AtomQuery query);

/**
 * @return the atomic symbol for the given Element enum
 */
SKETCHER_API QString get_atomic_symbol_from_element(Element element);

/**
 * @return a qstring that contains a list of all atomic symbols in the input,
 * separated by separator
 */
SKETCHER_API QString join_all_atomic_symbols(std::vector<Element> elements,
                                             QString separator = ",");

} // namespace sketcher
} // namespace schrodinger
