#pragma once

#include <string>
#include <tuple>
#include <unordered_map>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/QueryBond.h>
#include <rdkit/GraphMol/QueryOps.h>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

enum class BondTopology { IN_RING, NOT_IN_RING, EITHER };

// map of bond tool to {bond type, bond stereochemistry, whether the bond is
// flippable (i.e. does the bond change meaningfully if we swap the start and
// end atoms), icon for the cursor hint}
const std::unordered_map<
    BondTool,
    std::tuple<RDKit::Bond::BondType, RDKit::Bond::BondDir, bool, QString>>
    BOND_TOOL_BOND_MAP = {
        {BondTool::SINGLE,
         {RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::NONE, false,
          ":/icons/bond_single.svg"}},
        {BondTool::DOUBLE,
         {RDKit::Bond::BondType::DOUBLE, RDKit::Bond::BondDir::NONE, false,
          ":/icons/bond_double.svg"}},
        {BondTool::TRIPLE,
         {RDKit::Bond::BondType::TRIPLE, RDKit::Bond::BondDir::NONE, false,
          ":/icons/bond_triple.svg"}},
        {BondTool::COORDINATE,
         {RDKit::Bond::BondType::DATIVE, RDKit::Bond::BondDir::NONE, true,
          ":/icons/bond_coordinate.svg"}},
        {BondTool::ZERO,
         {RDKit::Bond::BondType::ZERO, RDKit::Bond::BondDir::NONE, false,
          ":/icons/bond_zero.svg"}},
        {BondTool::SINGLE_UP,
         {RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::BEGINWEDGE, true,
          ":/icons/bond_up.svg"}},
        {BondTool::SINGLE_DOWN,
         {RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::BEGINDASH, true,
          ":/icons/bond_down.svg"}},
        {BondTool::AROMATIC,
         {RDKit::Bond::BondType::AROMATIC, RDKit::Bond::BondDir::NONE, false,
          ":/icons/bond_aromatic.svg"}},
        {BondTool::SINGLE_EITHER,
         {RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::UNKNOWN, true,
          ":/icons/bond_wiggly.svg"}},
        {BondTool::DOUBLE_EITHER,
         {RDKit::Bond::BondType::DOUBLE, RDKit::Bond::BondDir::EITHERDOUBLE,
          true, ":/icons/bond_crossed.svg"}},
};

const std::unordered_map<AtomQuery,
                         std::function<RDKit::QueryAtom::QUERYATOM_QUERY*()>>
    ATOM_TOOL_QUERY_MAP = {
        {AtomQuery::A, RDKit::makeAAtomQuery},
        {AtomQuery::AH, RDKit::makeAHAtomQuery},
        {AtomQuery::Q, RDKit::makeQAtomQuery},
        {AtomQuery::QH, RDKit::makeQHAtomQuery},
        {AtomQuery::M, RDKit::makeMAtomQuery},
        {AtomQuery::MH, RDKit::makeMHAtomQuery},
        {AtomQuery::X, RDKit::makeXAtomQuery},
        {AtomQuery::XH, RDKit::makeXHAtomQuery},
};

const std::unordered_map<
    BondTool, std::pair<std::function<RDKit::QueryBond::QUERYBOND_QUERY*()>,
                        RDKit::Bond::BondType>>
    BOND_TOOL_QUERY_MAP = {
        {BondTool::SINGLE_OR_DOUBLE,
         {RDKit::makeSingleOrDoubleBondQuery, RDKit::Bond::BondType::SINGLE}},
        {BondTool::SINGLE_OR_AROMATIC,
         {RDKit::makeSingleOrAromaticBondQuery, RDKit::Bond::BondType::SINGLE}},
        {BondTool::DOUBLE_OR_AROMATIC,
         {RDKit::makeDoubleOrAromaticBondQuery, RDKit::Bond::BondType::DOUBLE}},
        {BondTool::ANY,
         {RDKit::makeBondNullQuery, RDKit::Bond::BondType::SINGLE}},
};

/**
 * For a given bond, returns a pair of
 *   - the type of bond we should draw this as, which will be SINGLE for all
 *     bonds with a query label
 *   - the appropriate query label for the bond
 * @param bond The bond to examine.  Must be non-null.
 */
SKETCHER_API std::pair<RDKit::Bond::BondType, std::string>
get_bond_type_and_query_label(const RDKit::Bond* const bond);

/**
 * Get the label string to display for the given bond query.  Returns "Query" if
 * the query cannot be parsed.
 * @param query The query to parse.  Must be non-null.
 */
SKETCHER_API std::string
get_label_for_bond_query(const RDKit::Bond::QUERYBOND_QUERY* const query);

/**
 * @return whether there's at least one implicit H on any of the specified
 * atoms.
 * @param atoms The atoms to consider.
 * @note This method is used to determine whether the set needs to have
 * explicit Hs added or removed.
 */
SKETCHER_API bool
has_any_implicit_Hs(const std::unordered_set<const RDKit::Atom*>& atoms);

/**
 * @return the topology of the bond.
 * @note This assumes that if the topology is set the
 * bond is a composite and between topology and something else (usually bond
 * order). This is enough to recognize topology queries created with the
 * sketcher, but may not be enough for all cases.
 */
SKETCHER_API BondTopology get_bond_topology(const RDKit::Bond* const bond);

/**
 * Set the topology of the bond.
 * @param bond The bond to modify.  Must be a query bond already.
 * @param topology The new topology for the bond. If the bond already has a
 * topology it will be updated to the new value
 */
SKETCHER_API void set_bond_topology(RDKit::QueryBond* const bond,
                                    BondTopology topology);

/**
 * Creates a new bond based on an existing query bond, removing the topology
 * info, assuming that the query is a BondAnd query between topology and
 * something else. If the existing bond does not have a query or the query is
 * not a "BondAnd" query or has more than two children, the new bond will be a
 * copy of the existing bond. If the remaining  query is a "BondOrder" query, a
 * new bond of the corresponding order will be returned.
 *
 * @param existing_bond A pointer to the existing RDKit::QueryBond object.
 * @return A shared pointer to the new RDKit::Bond object.
 */
SKETCHER_API std::shared_ptr<RDKit::Bond>
make_new_bond_without_topology(const RDKit::QueryBond* existing_bond);

} // namespace sketcher
} // namespace schrodinger
