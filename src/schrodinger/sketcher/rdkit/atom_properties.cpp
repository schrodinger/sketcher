#include "schrodinger/sketcher/rdkit/atom_properties.h"

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/QueryAtom.h>

namespace schrodinger
{
namespace sketcher
{

bool AbstractAtomProperties::operator==(
    const AbstractAtomProperties& other) const
{
    if (isQuery() != other.isQuery() || this->element != other.element ||
        this->isotope != other.isotope ||
        this->enhanced_stereo_type != other.enhanced_stereo_type ||
        this->enhanced_stereo_group_id != other.enhanced_stereo_group_id) {
        return false;
    } else if (isQuery()) {
        auto* this_query = static_cast<const AtomQueryProperties*>(this);
        auto* other_query = static_cast<const AtomQueryProperties*>(&other);
        return this_query->charge == other_query->charge &&
               this_query->unpaired_electrons ==
                   other_query->unpaired_electrons &&
               this_query->query_type == other_query->query_type &&
               this_query->allowed_list == other_query->allowed_list &&
               this_query->wildcard == other_query->wildcard &&
               this_query->r_group == other_query->r_group &&
               this_query->total_h == other_query->total_h &&
               this_query->num_connections == other_query->num_connections &&
               this_query->aromaticity == other_query->aromaticity &&
               this_query->ring_count_type == other_query->ring_count_type &&
               this_query->ring_count_exact_val ==
                   other_query->ring_count_exact_val &&
               this_query->ring_bond_count_type ==
                   other_query->ring_bond_count_type &&
               this_query->ring_bond_count_exact_val ==
                   other_query->ring_bond_count_exact_val &&
               this_query->smallest_ring_size ==
                   other_query->smallest_ring_size &&
               this_query->smarts_query == other_query->smarts_query;
    } else { // not a query
        auto* this_atom = static_cast<const AtomProperties*>(this);
        auto* other_atom = static_cast<const AtomProperties*>(&other);
        return this_atom->charge == other_atom->charge &&
               this_atom->unpaired_electrons == other_atom->unpaired_electrons;
    }
}

std::shared_ptr<AbstractAtomProperties>
read_properties_from_atom(const RDKit::Atom* const atom)
{
    // TODO
    return std::make_shared<AtomProperties>(AtomProperties());
}

RDKit::Atom create_atom_with_properties(
    const std::shared_ptr<AbstractAtomProperties> properties)
{
    // TODO
    if (properties->isQuery()) {
        // auto* query_props =
        //     static_cast<const AtomQueryProperties*>(properties.get());
        return RDKit::QueryAtom(0);
    } else {
        // auto* atom_props = static_cast<const
        // AtomProperties*>(properties.get());
        return RDKit::Atom(6);
    }
}

} // namespace sketcher
} // namespace schrodinger
