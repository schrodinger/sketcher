#include "schrodinger/sketcher/rdkit/monomeric.h"

#include <functional>

#include <QColor>
#include <QGraphicsItem>
#include <QPointF>

#include <fmt/core.h>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{
const std::string PEPTIDE_POLYMER_PREFIX = "PEPTIDE";
// According to HELM, DNA is a subtype of RNA, so DNA also uses the RNA prefix
const std::string NUCLEOTIDE_POLYMER_PREFIX = "RNA";

// "Pretty" names for attachment points that are normally represented as "R#"
// Note that the primes use apostrophes instead of a Unicode prime to avoid
// issues with C++′s handling of Unicode.
const std::unordered_map<MonomerType, std::vector<std::string>>
    NUMBERED_AP_NAMES = {
        {MonomerType::PEPTIDE, {"N", "C", "X"}},
        {MonomerType::NA_BASE, {"N1/9"}},
        {MonomerType::NA_SUGAR, {"3'", "5'", "1'"}},
};

// standard attachment point names that don't follow the "R#" naming scheme
const std::unordered_map<MonomerType, std::vector<std::string>>
    EXPECTED_AP_CUSTOM_NAMES = {{MonomerType::NA_BASE, {"pair"}}};

const int INVALID_ATTACHMENT_POINT_SPEC = -1;
const int ATTACHMENT_POINT_WITH_CUSTOM_NAME = -2;

} // namespace

MonomerType get_monomer_type(const RDKit::Atom* atom)
{
    const auto* monomer_info = atom->getMonomerInfo();
    if (monomer_info == nullptr) {
        throw std::runtime_error("Atom has no monomer info");
    }
    const auto* res_info =
        dynamic_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info);
    if (res_info == nullptr) {
        return MonomerType::CHEM;
    }
    const auto& chain_id = res_info->getChainId();
    if (chain_id.starts_with(PEPTIDE_POLYMER_PREFIX)) {
        return MonomerType::PEPTIDE;
    } else if (chain_id.starts_with(NUCLEOTIDE_POLYMER_PREFIX)) {
        const auto& res_name = res_info->getResidueName();
        if (res_name.empty()) {
            return MonomerType::NA_BASE;
        }
        auto last_char = std::tolower(res_name.back());
        if (last_char == 'p') {
            return MonomerType::NA_PHOSPHATE;
        } else if (last_char == 'r') {
            return MonomerType::NA_SUGAR;
        }
        return MonomerType::NA_BASE;
    }
    return MonomerType::CHEM;
}

std::string get_monomer_res_name(const RDKit::Atom* const monomer)
{
    bool is_smiles = false;
    monomer->getPropIfPresent(SMILES_MONOMER, is_smiles);
    if (is_smiles) {
        return SMILES_PLACEHOLDER_TEXT;
    }
    const auto* monomer_info = monomer->getMonomerInfo();
    const auto* res_info =
        dynamic_cast<const RDKit::AtomPDBResidueInfo*>(monomer_info);
    if (res_info == nullptr) {
        return monomer_info->getName();
    }
    return res_info->getResidueName();
}

bool contains_two_monomer_linkages(const RDKit::Bond* bond)
{
    std::string linkage, custom_linkage;
    bond->getPropIfPresent(LINKAGE, linkage);
    bool custom_linkage_exists =
        bond->getPropIfPresent(CUSTOM_BOND, custom_linkage);
    return custom_linkage_exists && custom_linkage != linkage;
}

std::tuple<bool, bool, QColor, qreal, Qt::PenStyle>
get_connector_style(const RDKit::Bond* bond, const bool is_secondary_connection)
{
    const auto* start_atom = bond->getBeginAtom();
    auto start_res_name = get_monomer_res_name(start_atom);
    auto start_monomer_type = get_monomer_type(start_atom);

    const auto* end_atom = bond->getEndAtom();
    auto end_res_name = get_monomer_res_name(end_atom);
    auto end_monomer_type = get_monomer_type(end_atom);

    if (start_monomer_type == MonomerType::CHEM ||
        end_monomer_type == MonomerType::CHEM) {
        // chem connector
        return {false, false, CHEM_CONNECTOR_COLOR, CHEM_CONNECTOR_WIDTH,
                Qt::PenStyle::SolidLine};
    } else if (start_monomer_type == MonomerType::PEPTIDE ||
               end_monomer_type == MonomerType::PEPTIDE) {
        std::string attachment_points;
        std::string prop = is_secondary_connection ? CUSTOM_BOND : LINKAGE;
        bond->getPropIfPresent(prop, attachment_points);
        if (start_res_name.ends_with('C') && end_res_name.ends_with('C') &&
            attachment_points == "R3-R3") {
            return {true, true, DISULFIDE_CONNECTOR_COLOR,
                    DISULFIDE_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine};
        }
        bool start_is_branch = false;
        bool end_is_branch = false;
        start_atom->getPropIfPresent(BRANCH_MONOMER, start_is_branch);
        end_atom->getPropIfPresent(BRANCH_MONOMER, end_is_branch);
        if (start_is_branch || end_is_branch) {
            // branching connector
            return {start_is_branch, end_is_branch,
                    AA_BRANCHING_CONNECTOR_COLOR, AA_BRANCHING_CONNECTOR_WIDTH,
                    Qt::PenStyle::SolidLine};
        }
        // standard amino acid linear connector
        return {false, false, AA_LINEAR_CONNECTOR_COLOR,
                AA_LINEAR_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine};
    } else if (start_monomer_type == MonomerType::NA_BASE &&
               end_monomer_type == MonomerType::NA_BASE) {
        // nucleic acid base connector
        return {false, false, NA_BASE_CONNECTOR_COLOR, NA_BASE_CONNECTOR_WIDTH,
                Qt::PenStyle::DotLine};
    } else if (start_monomer_type != MonomerType::NA_BASE &&
               end_monomer_type != MonomerType::NA_BASE) {
        // nucleic acid backbone connector
        return {false, false, NA_BACKBONE_CONNECTOR_COLOR,
                NA_BACKBONE_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine};
    } else {
        // nucleic acid backbone to base connector
        return {false, false, NA_BACKBONE_TO_BASE_CONNECTOR_COLOR,
                NA_BACKBONE_TO_BASE_CONNECTOR_WIDTH, Qt::PenStyle::SolidLine};
    }
}

/**
 * Return true if coord is above (within a 45 degree cone of) other
 */
static bool is_coord_above_the_other(const QPointF& coord, const QPointF& other)
{
    return coord.y() < other.y() &&
           qFabs(coord.x() - other.x()) < qFabs(coord.y() - other.y());
}

qreal get_monomer_arrowhead_offset(const QGraphicsItem& monomer_item,
                                   const QPointF& bound_coords)
{
    auto offset = monomer_item.boundingRect().height() / 2 +
                  MONOMER_CONNECTOR_ARROWHEAD_RADIUS;
    if (is_coord_above_the_other(monomer_item.pos(), bound_coords)) {
        offset *= -1;
    }
    return offset;
}

/**
 * Return the name and number of the attachment point specified in the given
 * linkage string.
 * @param linkage A description of the linkage taken from the bond properties.
 * Formatted similar to "R2-R3".
 * @param is_begin_atom If true, we'll return the attachment point for the
 * bond's begin atom.  Otherwise, we'll return the attachment point for the
 * bond's end atom.
 * @return A pair of
 *   - The attachment point number. For attachment points named as R# (which is
 *     most attachment points), this will be the given number. For attachment
 *     points with a custom name (e.g. "pair"), this number will be
 *     ATTACHMENT_POINT_WITH_CUSTOM_NAME. If the attachment point specification
 *     could not be parsed (i.e. no dash was found), this number will be
 *     INVALID_ATTACHMENT_POINT_SPEC.
 *   - The attachment point name (e.g. "R1"), or empty string if the attachment
 *     point specification could not be parsed.
 */
static std::pair<int, std::string>
get_attachment_point_for_atom(std::string linkage, bool is_begin_atom)
{
    auto dash_pos = linkage.find("-");
    if (dash_pos == std::string::npos) {
        return {INVALID_ATTACHMENT_POINT_SPEC, ""};
    }
    std::string attachment_point_name = is_begin_atom
                                            ? linkage.substr(0, dash_pos)
                                            : linkage.substr(dash_pos + 1);
    if (attachment_point_name[0] != 'R') {
        return {ATTACHMENT_POINT_WITH_CUSTOM_NAME, attachment_point_name};
    }
    // remove the leading 'R' now that we've confirmed it exists
    attachment_point_name.erase(0, 1);
    int ap_num;
    try {
        ap_num = std::stoi(attachment_point_name);
    } catch (const std::logic_error&) {
        // it's not an integer
        ap_num = ATTACHMENT_POINT_WITH_CUSTOM_NAME;
    }
    return {ap_num, attachment_point_name};
}

/**
 * @overload Return the number of the attachment point for the bond on the
 * specified monomer. Note that if the bond specifies two linkages, this
 * overload will only return the attachment point of the primary linkage.
 * @return the attachment point number, or -1 if the attachment point is not
 * properly specified
 */
static std::pair<int, std::string>
get_attachment_point_for_atom(const RDKit::Atom* monomer,
                              const RDKit::Bond* bond)
{
    std::string linkage;
    if (bond->getPropIfPresent(LINKAGE, linkage)) {
        bool is_start_atom = bond->getBeginAtom() == monomer;
        return get_attachment_point_for_atom(linkage, is_start_atom);
    }
    return {INVALID_ATTACHMENT_POINT_SPEC, ""};
}

/**
 * @return a list of information about all bound attachment points of the
 * specified monomer. The following information is given for each attachment
 * points:
 *   - The number of the attachment point, e.g., 1 for "R1". Will be
 *     ATTACHMENT_POINT_WITH_CUSTOM_NAME if the attachment point uses a custom
 *     name (e.g. "pair")
 *   - The name of the attachment point, e.g. "R1"
 *   - The other monomer involved in the connection
 *   - Whether this attachment point is involved in the secondary connection of
 *     the bond. Secondary connections occur when a single RDKit::Bond
 *     represents multiple connections, e.g. two neighboring cysteines that are
 *     disulfide bonded to each other.
 */
static std::vector<std::tuple<int, std::string, const RDKit::Atom*, bool>>
get_bound_attachment_points(const RDKit::Atom* monomer)
{
    const auto& mol = monomer->getOwningMol();
    std::unordered_set<int> bound_ap_nums;
    std::unordered_set<std::string> bound_ap_custom_names;
    std::vector<std::tuple<int, std::string, const RDKit::Atom*, bool>>
        bound_aps;

    auto record_linkage = [&bound_aps, &bound_ap_nums, &bound_ap_custom_names,
                           monomer](const RDKit::Bond* bond,
                                    const bool is_start_atom,
                                    const std::string& prop_name) {
        std::string linkage;
        if (bond->getPropIfPresent(prop_name, linkage)) {
            const auto& [ap_value, ap_name] =
                get_attachment_point_for_atom(linkage, is_start_atom);
            if (ap_value > 0 && !bound_ap_nums.contains(ap_value)) {
                bound_ap_nums.insert(ap_value);
            } else if (ap_value == ATTACHMENT_POINT_WITH_CUSTOM_NAME &&
                       !bound_ap_custom_names.contains(ap_name)) {
                bound_ap_custom_names.insert(ap_name);
            } else {
                return;
            }
            bound_aps.push_back({ap_value, ap_name, bond->getOtherAtom(monomer),
                                 prop_name == CUSTOM_BOND});
        }
    };

    for (auto* bond : mol.atomBonds(monomer)) {
        bool is_start_atom = monomer == bond->getBeginAtom();
        record_linkage(bond, is_start_atom, LINKAGE);
        record_linkage(bond, is_start_atom, CUSTOM_BOND);
    }
    return bound_aps;
}

/**
 * @return Sets of all attachment points on the specified monomer that are
 * currently unbound. Two sets are returned:
 *   - Numbered attachment points. Attachment points are specified using
 *     integers, e.g. 1 for "R1".
 *   - Attachment points with custom names (e.g. "pair")
 *
 * Note that we don't have a good way to determine how many attachment points a
 * CHEM monomer should have, so we assume that it has one additional attachment
 * point beyond the highest numbered bound attachment point.
 */
static std::pair<std::vector<int>, std::vector<std::string>>
get_available_attachment_points(const RDKit::Atom* monomer)
{
    auto bound_aps = get_bound_attachment_points(monomer);
    auto monomer_type = get_monomer_type(monomer);
    std::unordered_set<int> bound_ap_nums;
    std::unordered_set<std::string> bound_aps_with_custom_names;
    for (const auto& [ap_num, orig_ap_name, atom, is_secondary_connection] :
         bound_aps) {
        if (ap_num > 0) {
            bound_ap_nums.insert(ap_num);
        } else if (ap_num == ATTACHMENT_POINT_WITH_CUSTOM_NAME) {
            bound_aps_with_custom_names.insert(orig_ap_name);
        }
    }
    int num_numbered_aps = -1;
    if (NUMBERED_AP_NAMES.contains(monomer_type)) {
        num_numbered_aps = NUMBERED_AP_NAMES.at(monomer_type).size();
    } else if (monomer_type == MonomerType::NA_PHOSPHATE) {
        num_numbered_aps = 2;
    } else {
        // a CHEM monomer
        num_numbered_aps =
            *std::max_element(bound_ap_nums.begin(), bound_ap_nums.end());
        num_numbered_aps += 1;
    }
    std::vector<int> available_numbered_aps;
    for (int ap = 1; ap <= num_numbered_aps; ++ap) {
        if (!bound_ap_nums.contains(ap)) {
            available_numbered_aps.push_back(ap);
        }
    }
    std::vector<std::string> available_aps_with_custom_names;
    if (EXPECTED_AP_CUSTOM_NAMES.contains(monomer_type)) {
        for (const auto& ap_name : EXPECTED_AP_CUSTOM_NAMES.at(monomer_type)) {
            if (!bound_aps_with_custom_names.contains(ap_name)) {
                available_aps_with_custom_names.push_back(ap_name);
            }
        }
    }
    return {available_numbered_aps, available_aps_with_custom_names};
}

/**
 * Convert an numbered attachment point (i.e. one named "R#") number to a
 * "pretty" name
 * @param ap_num The attachment point number to convert. Should be > 0, as
 * attachment points are 1-indexed.
 * @param all_names A list of "pretty" names for attachment points, starting
 * with R1.
 * @return If all_names contains a "pretty" name for ap_num, then that name will
 * be returned. Otherwise "R<ap_num>" will be returned.
 */
static std::string ap_num_to_name(const int ap_num,
                                  const std::vector<std::string>& all_names)
{
    if (0 < ap_num && static_cast<unsigned int>(ap_num) <= all_names.size()) {
        return all_names[ap_num - 1];
    }
    return fmt::format("R{}", ap_num);
}

/**
 * If the given phosphate monomer is bound to exactly one sugar (or if it's
 * bound to a chain of phosphates and that chain of phosphates is bound to a
 * sugar, e.g. ATP), return the "pretty" name of the sugar's attachment point
 * (e.g. "3'", not "R1"). Otherwise, return en empty string.
 */
static std::string
get_attachment_point_name_of_bound_sugar(const RDKit::Atom* phosphate)
{
    const auto& mol = phosphate->getOwningMol();
    if (mol.getAtomDegree(phosphate) != 1) {
        return "";
    }
    auto prev_neighbor = phosphate;
    auto cur_neighbor = *mol.atomNeighbors(phosphate).begin();
    // if there's a chain of phosphates, continue along it until we reach the
    // sugar
    while (mol.getAtomDegree(cur_neighbor) == 2 &&
           get_monomer_type(cur_neighbor) == MonomerType::NA_PHOSPHATE) {
        for (auto possible_next_neighbor : mol.atomNeighbors(cur_neighbor)) {
            if (possible_next_neighbor != prev_neighbor) {
                prev_neighbor = cur_neighbor;
                cur_neighbor = possible_next_neighbor;
                break;
            }
        }
    }
    if (get_monomer_type(cur_neighbor) == MonomerType::NA_SUGAR) {
        auto bond_to_sugar = mol.getBondBetweenAtoms(prev_neighbor->getIdx(),
                                                     cur_neighbor->getIdx());
        const auto& [sugar_ap_num, ap_name] =
            get_attachment_point_for_atom(cur_neighbor, bond_to_sugar);
        // the phosphate should be bound to either the 3' (R1) or 5' (R2). If
        // it's bound to something else, ignore it since something's gone wrong.
        if ((sugar_ap_num == 1 || sugar_ap_num == 2)) {
            return ap_num_to_name(sugar_ap_num,
                                  NUMBERED_AP_NAMES.at(MonomerType::NA_SUGAR));
        }
    }
    return "";
}

/**
 * @return a list of all "pretty" attachment point names (e.g. "N" instead of
 * "R1" for amino acids) for numbered attachment points (i.e. any attachment
 * points named "R#") of the given monomer.  This function does not account for
 * of whether those attachment points are bound or available.
 *
 * Note that CHEM monomers don't have special names, so we return an empty list
 * (which will cause ap_num_to_name() to return R1, R2, etc).
 *
 * Also note that phosphate attachment point names reflect the attachment point
 * of the bound sugar, as the sites themselves are chemically identical. As a
 * result, these names are only meaningful when exactly one sugar is bound.
 * Because of this, we return blank attachment point names (i.e. empty strings,
 * *not* an empty list) unless there is exactly one attachment point bound to a
 * sugar.
 */
static std::vector<std::string>
get_all_numbered_attachment_point_names(const RDKit::Atom* monomer)
{
    auto monomer_type = get_monomer_type(monomer);

    if (NUMBERED_AP_NAMES.contains(monomer_type)) {
        return NUMBERED_AP_NAMES.at(monomer_type);
    } else if (monomer_type == MonomerType::NA_PHOSPHATE) {
        std::vector<std::string> phos_ap_names = {"", ""};
        auto sugar_ap_name = get_attachment_point_name_of_bound_sugar(monomer);
        if (!sugar_ap_name.empty()) {
            const auto& mol = monomer->getOwningMol();
            // the phosphate must have exactly one bond; otherwise,
            // sugar_ap_name would be empty
            const RDKit::Bond* phos_bond = *mol.atomBonds(monomer).begin();
            const auto& [bound_phos_ap_num, ap_name] =
                get_attachment_point_for_atom(monomer, phos_bond);
            if (bound_phos_ap_num == 1 || bound_phos_ap_num == 2) {
                // phosphates should only have two attachment points
                int unbound_phos_ap_name_idx = bound_phos_ap_num == 1 ? 1 : 0;
                phos_ap_names[unbound_phos_ap_name_idx] = sugar_ap_name;
            }
        }
        return phos_ap_names;
    } else {
        // for CHEM monomers, we return an empty list, meaning that the
        // attachment points will be named R1, R2, etc
        return {};
    }
}

std::vector<std::tuple<std::string, const RDKit::Atom*, bool>>
get_bound_attachment_point_names_and_atoms(const RDKit::Atom* monomer)
{
    auto bound_aps = get_bound_attachment_points(monomer);
    auto all_names = get_all_numbered_attachment_point_names(monomer);
    std::vector<std::tuple<std::string, const RDKit::Atom*, bool>>
        bound_ap_info;
    std::transform(
        bound_aps.begin(), bound_aps.end(), std::back_inserter(bound_ap_info),
        [&all_names](auto ap_info) {
            const auto& [ap_num, orig_ap_name, atom, is_secondary_connection] =
                ap_info;
            auto pretty_ap_name = ap_num == ATTACHMENT_POINT_WITH_CUSTOM_NAME
                                      ? orig_ap_name
                                      : ap_num_to_name(ap_num, all_names);
            return std::make_tuple(pretty_ap_name, atom,
                                   is_secondary_connection);
        });
    return bound_ap_info;
}

std::vector<std::string>
get_available_attachment_point_names(const RDKit::Atom* monomer)
{
    auto [available_numbered_aps, available_custom_named_aps] =
        get_available_attachment_points(monomer);
    auto all_names = get_all_numbered_attachment_point_names(monomer);
    std::vector<std::string> available_names;
    std::transform(available_numbered_aps.begin(), available_numbered_aps.end(),
                   std::back_inserter(available_names),
                   std::bind(ap_num_to_name, std::placeholders::_1, all_names));
    std::copy(available_custom_named_aps.begin(),
              available_custom_named_aps.end(),
              std::back_inserter(available_names));
    return available_names;
}

std::string
get_attachment_point_name_for_atom(const RDKit::Atom* monomer,
                                   const RDKit::Bond* connector,
                                   const bool is_secondary_connection)
{
    auto all_names = get_all_numbered_attachment_point_names(monomer);
    std::string prop_name = is_secondary_connection ? CUSTOM_BOND : LINKAGE;
    std::string linkage;
    if (connector->getPropIfPresent(prop_name, linkage)) {
        bool is_start_atom = connector->getBeginAtom() == monomer;
        const auto& [ap_num, orig_ap_name] =
            get_attachment_point_for_atom(linkage, is_start_atom);
        if (ap_num > 0) {
            return ap_num_to_name(ap_num, all_names);
        } else if (ap_num == ATTACHMENT_POINT_WITH_CUSTOM_NAME) {
            return orig_ap_name;
        }
    }
    return "";
}

} // namespace sketcher
} // namespace schrodinger
