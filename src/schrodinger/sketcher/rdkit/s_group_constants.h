#pragma once
#include <string>

#include <boost/assign.hpp>
#include <boost/bimap.hpp>

namespace schrodinger
{
namespace sketcher
{

/**
 * Subgroup types for bracket subgroups
 */
enum class SubgroupType {
    SRU_POLYMER,
    COPOLYMER,
    OTHER,
};

/**
 * Repeat patterns for bracket subgroups
 */
enum class RepeatPattern {
    EITHER_UNKNOWN,
    HEAD_TO_HEAD,
    HEAD_TO_TAIL,
    UNSPECIFIED,
};

using RepeatPatternStringBimapType = boost::bimap<RepeatPattern, std::string>;
const RepeatPatternStringBimapType REPEATPATTERN_TO_RDKITSTRING_BIMAP =
    boost::assign::list_of<RepeatPatternStringBimapType::relation>(
        RepeatPattern::HEAD_TO_TAIL, "HT")(RepeatPattern::HEAD_TO_HEAD, "HH")(
        RepeatPattern::EITHER_UNKNOWN, "EU");

// Currently, we only support 2 out of 15 SGroup types, so for now it makes
// sense to have a OTHER SubgroupType, which uses the GENeric SGroup type.
using SubgroupTypeStringBimapType = boost::bimap<SubgroupType, std::string>;
const SubgroupTypeStringBimapType SUBGROUPTYPE_TO_RDKITSTRING_BIMAP =
    boost::assign::list_of<SubgroupTypeStringBimapType::relation>(
        SubgroupType::SRU_POLYMER, "SRU")(SubgroupType::COPOLYMER,
                                          "COP")(SubgroupType::OTHER, "GEN");

} // namespace sketcher
} // namespace schrodinger
