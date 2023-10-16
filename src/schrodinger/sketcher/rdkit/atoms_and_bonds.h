
#pragma once

#include <string>
#include <tuple>
#include <unordered_map>

#include <GraphMol/Atom.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <GraphMol/QueryOps.h>

#include "schrodinger/sketcher/model/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

// map of bond tool to {bond type, bond stereochemistry, whether the bond is
// flippable (i.e. does the bond change meaningfully if we swap the start
// and end atoms)}
const std::unordered_map<
    BondTool, std::tuple<RDKit::Bond::BondType, RDKit::Bond::BondDir, bool>>
    BOND_TOOL_BOND_MAP = {
        {BondTool::SINGLE,
         {RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::NONE, false}},
        {BondTool::DOUBLE,
         {RDKit::Bond::BondType::DOUBLE, RDKit::Bond::BondDir::NONE, false}},
        {BondTool::TRIPLE,
         {RDKit::Bond::BondType::TRIPLE, RDKit::Bond::BondDir::NONE, false}},
        {BondTool::COORDINATE,
         {RDKit::Bond::BondType::DATIVE, RDKit::Bond::BondDir::NONE, true}},
        {BondTool::ZERO,
         {RDKit::Bond::BondType::ZERO, RDKit::Bond::BondDir::NONE, false}},
        {BondTool::SINGLE_UP,
         {RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::BEGINWEDGE,
          true}},
        {BondTool::SINGLE_DOWN,
         {RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::BEGINDASH,
          true}},
        {BondTool::AROMATIC,
         {RDKit::Bond::BondType::AROMATIC, RDKit::Bond::BondDir::NONE, false}},
        {BondTool::SINGLE_EITHER,
         {RDKit::Bond::BondType::SINGLE, RDKit::Bond::BondDir::UNKNOWN, true}},
        {BondTool::DOUBLE_EITHER,
         {RDKit::Bond::BondType::DOUBLE, RDKit::Bond::BondDir::EITHERDOUBLE,
          true}},
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

// map of bond tool to {the label text, a function that returns the appropriate
// query object)}
const std::unordered_map<
    BondTool,
    std::pair<std::string, std::function<RDKit::QueryBond::QUERYBOND_QUERY*()>>>
    BOND_TOOL_QUERY_MAP = {{BondTool::SINGLE_OR_DOUBLE,
                            {"S/D", RDKit::makeSingleOrDoubleBondQuery}},
                           {BondTool::SINGLE_OR_AROMATIC,
                            {"S/A", RDKit::makeSingleOrAromaticBondQuery}},
                           {BondTool::DOUBLE_OR_AROMATIC,
                            {"D/A", RDKit::makeDoubleOrAromaticBondQuery}},
                           {BondTool::ANY, {"Any", RDKit::makeBondNullQuery}}};

} // namespace sketcher
} // namespace schrodinger
