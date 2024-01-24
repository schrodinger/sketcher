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

const std::unordered_map<BondTool,
                         std::function<RDKit::QueryBond::QUERYBOND_QUERY*()>>
    BOND_TOOL_QUERY_MAP = {
        {BondTool::SINGLE_OR_DOUBLE, RDKit::makeSingleOrDoubleBondQuery},
        {BondTool::SINGLE_OR_AROMATIC, RDKit::makeSingleOrAromaticBondQuery},
        {BondTool::DOUBLE_OR_AROMATIC, RDKit::makeDoubleOrAromaticBondQuery},
        {BondTool::ANY, RDKit::makeBondNullQuery}};

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

} // namespace sketcher
} // namespace schrodinger
