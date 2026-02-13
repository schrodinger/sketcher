#include "schrodinger/sketcher/rdkit/monomeric.h"

#include <functional>
#include <optional>

#include <boost/range/join.hpp>

#include <QGraphicsItem>
#include <QPointF>

#include <fmt/core.h>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/molviewer/monomer_constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

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
    NUMBERED_AP_NAMES_BY_MONOMER_TYPE = {
        {MonomerType::PEPTIDE, {"N", "C", "X"}},
        {MonomerType::NA_BASE, {"N1/9"}},
        {MonomerType::NA_SUGAR, {"3'", "5'", "1'"}},
};

// standard attachment point names that don't follow the "R#" naming scheme
const std::unordered_map<MonomerType, std::vector<std::string>>
    EXPECTED_AP_CUSTOM_NAMES = {{MonomerType::NA_BASE, {"pair"}}};

const int INVALID_ATTACHMENT_POINT_SPEC = -2;

const std::unordered_map<Direction, Direction> OPPOSITE_CARDINAL_DIRECTION = {
    {Direction::N, Direction::S},
    {Direction::S, Direction::N},
    {Direction::E, Direction::W},
    {Direction::W, Direction::E}};
const std::unordered_map<Direction, Direction> CLOCKWISE_CARDINAL_DIRECTION = {
    {Direction::N, Direction::E},
    {Direction::S, Direction::W},
    {Direction::E, Direction::S},
    {Direction::W, Direction::N}};
const std::unordered_map<Direction, std::vector<Direction>>
    PERPENDICULAR_CARDINAL_DIRECTIONS = {
        {Direction::N, {Direction::W, Direction::E}},
        {Direction::S, {Direction::W, Direction::E}},
        {Direction::E, {Direction::N, Direction::S}},
        {Direction::W, {Direction::N, Direction::S}}};

class NoAvailableDirectionsException : public std::exception
{
};

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

ConnectorType get_connector_type(const RDKit::Bond* bond,
                                 const bool is_secondary_connection)
{
    const auto* start_atom = bond->getBeginAtom();
    auto start_res_name = get_monomer_res_name(start_atom);
    auto start_monomer_type = get_monomer_type(start_atom);

    const auto* end_atom = bond->getEndAtom();
    auto end_res_name = get_monomer_res_name(end_atom);
    auto end_monomer_type = get_monomer_type(end_atom);

    if (start_monomer_type == MonomerType::CHEM ||
        end_monomer_type == MonomerType::CHEM) {
        return ConnectorType::CHEM;
    } else if (start_monomer_type == MonomerType::PEPTIDE ||
               end_monomer_type == MonomerType::PEPTIDE) {
        std::string attachment_points;
        std::string prop = is_secondary_connection ? CUSTOM_BOND : LINKAGE;
        bond->getPropIfPresent(prop, attachment_points);
        if (attachment_points == "R3-R3") {
            if (start_res_name.ends_with('C') && end_res_name.ends_with('C')) {
                return ConnectorType::PEPTIDE_DISULFIDE;
            } else {
                return ConnectorType::PEPTIDE_SIDE_CHAIN;
            }
        }
        bool start_is_branch = false;
        bool end_is_branch = false;
        start_atom->getPropIfPresent(BRANCH_MONOMER, start_is_branch);
        end_atom->getPropIfPresent(BRANCH_MONOMER, end_is_branch);
        if (start_is_branch || end_is_branch) {
            return ConnectorType::PEPTIDE_BRANCHING;
        }
        return ConnectorType::PEPTIDE_LINEAR;
    } else if (start_monomer_type == MonomerType::NA_BASE &&
               end_monomer_type == MonomerType::NA_BASE) {
        return ConnectorType::NA_BASE;
    } else if (start_monomer_type != MonomerType::NA_BASE &&
               end_monomer_type != MonomerType::NA_BASE) {
        return ConnectorType::NA_BACKBONE;
    } else {
        return ConnectorType::NA_BACKBONE_TO_BASE;
    }
}

std::pair<bool, bool>
does_connector_have_arrowheads(const RDKit::Bond* bond,
                               const bool is_secondary_connection)
{
    auto connector_type = get_connector_type(bond, is_secondary_connection);
    return does_connector_have_arrowheads(bond, connector_type);
}

std::pair<bool, bool>
does_connector_have_arrowheads(const RDKit::Bond* bond,
                               const ConnectorType connector_type)
{
    switch (connector_type) {
        case ConnectorType::CHEM:
        case ConnectorType::PEPTIDE_LINEAR:
        case ConnectorType::NA_BASE:
        case ConnectorType::NA_BACKBONE:
        case ConnectorType::NA_BACKBONE_TO_BASE:
            return {false, false};

        case ConnectorType::PEPTIDE_DISULFIDE:
        case ConnectorType::PEPTIDE_SIDE_CHAIN:
            return {true, true};

        case ConnectorType::PEPTIDE_BRANCHING: {
            const auto* start_atom = bond->getBeginAtom();
            const auto* end_atom = bond->getEndAtom();
            bool start_is_branch = false;
            bool end_is_branch = false;
            start_atom->getPropIfPresent(BRANCH_MONOMER, start_is_branch);
            end_atom->getPropIfPresent(BRANCH_MONOMER, end_is_branch);
            return {start_is_branch, end_is_branch};
        }
        default:
            return {false, false};
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

static Direction cardinal_direction_for_point(const RDGeom::Point3D& point)
{
    if (std::fabs(point.x) >= std::fabs(point.y)) {
        return point.x > 0 ? Direction::E : Direction::W;
    } else {
        return point.y > 0 ? Direction::N : Direction::S;
    }
}

static Direction get_bound_attachment_point_cardinal_direction(
    const RDKit::Atom* const monomer, const RDKit::Atom* const bound_monomer,
    const bool is_secondary_connection)
{
    auto& mol = monomer->getOwningMol();
    auto& conf = mol.getConformer();
    auto monomer_coords = conf.getAtomPos(monomer->getIdx());
    auto bound_monomer_coords = conf.getAtomPos(bound_monomer->getIdx());

    auto* bond =
        mol.getBondBetweenAtoms(monomer->getIdx(), bound_monomer->getIdx());
    auto [is_arrowhead_at_bond_beginning, is_arrowhead_at_bond_end] =
        does_connector_have_arrowheads(bond, is_secondary_connection);
    bool is_start_atom = bond->getBeginAtom() == monomer;
    bool is_arrowhead = is_start_atom ? is_arrowhead_at_bond_beginning
                                      : is_arrowhead_at_bond_end;
    if (is_arrowhead) {
        bool is_above = is_coord_above_the_other(
            to_scene_xy(monomer_coords), to_scene_xy(bound_monomer_coords));
        return is_above ? Direction::S : Direction::N;
    } else {
        auto relative_pos = bound_monomer_coords - monomer_coords;
        return cardinal_direction_for_point(relative_pos);
    }
}

/**
 * @return a list of information about all bound attachment points of the
 * specified monomer. Note that directions will be determined using cardinal
 * directions only, not diagonals.
 */
static std::vector<BoundAttachmentPoint>
get_bound_attachment_points(const RDKit::Atom* monomer)
{
    const auto& mol = monomer->getOwningMol();
    std::unordered_set<int> bound_ap_nums;
    std::unordered_set<std::string> bound_ap_custom_names;
    std::vector<BoundAttachmentPoint> bound_aps;

    auto record_linkage = [&bound_aps, &bound_ap_nums, &bound_ap_custom_names,
                           monomer](const RDKit::Bond* bond,
                                    const bool is_start_atom,
                                    const std::string& prop_name) {
        std::string linkage;
        if (bond->getPropIfPresent(prop_name, linkage)) {
            const auto& [ap_num, ap_name] =
                get_attachment_point_for_atom(linkage, is_start_atom);
            if (ap_num > 0 && !bound_ap_nums.contains(ap_num)) {
                bound_ap_nums.insert(ap_num);
            } else if (ap_num == ATTACHMENT_POINT_WITH_CUSTOM_NAME &&
                       !bound_ap_custom_names.contains(ap_name)) {
                bound_ap_custom_names.insert(ap_name);
            } else {
                return;
            }
            auto bound_monomer = bond->getOtherAtom(monomer);
            bool is_secondary_connection = prop_name == CUSTOM_BOND;
            auto dir = get_bound_attachment_point_cardinal_direction(
                monomer, bound_monomer, is_secondary_connection);
            bound_aps.push_back(BoundAttachmentPoint(
                ap_name, ap_num, bound_monomer, is_secondary_connection, dir));
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
 * @return the direction associated with the specified numbered attachment point
 */
static std::optional<Direction> fetch_assigned_direction_for_attachment_point(
    const int ap_num, const std::vector<BoundAttachmentPoint>& bound_aps,
    const std::vector<UnboundAttachmentPoint>& unbound_aps)
{
    for (auto cur_ap : bound_aps) {
        if (cur_ap.num == ap_num) {
            return cur_ap.direction;
        }
    }
    for (auto cur_ap : unbound_aps) {
        if (cur_ap.num == ap_num) {
            return cur_ap.direction;
        }
    }
    return std::nullopt;
}

/**
 * Determine which direction the specified unbound attachment point should be
 * drawn in
 * @param ap_num The attachment point number of the attachment point to
 * determine the direction of, e.g. 3 for "R3".  Should be
 * ATTACHMENT_POINT_WITH_CUSTOM_NAME if the attachment point doesn't follow the
 * "R#" naming scheme, such as the "pair" attachment point of nucleic acid
 * bases.
 * @param ap_name The name of the attachment point to determine the direction
 * of, e.g. "R3" or "pair"
 * @param monomer_type The type of monomer that this attachment point is on
 * @param bound_aps A list of all bound attachment point for this monomer with
 * their assigned directions
 * @param unbound_aps A list of all unbound attachment points for this monomer
 * that have already had their direction assigned. For a numbered attachment
 * point, this function assumes this list includes all attachment points with
 * lower numbers (i.e. directions should be assigned in numerical order). For an
 * attachment point with a custom name (e.g. "pair"), this function assumes that
 * this list includes all numbered attachment points.
 * @param occupied_directions A set of all assigned directions found in
 * bound_aps and unbound_aps.
 * @return the direction to assign to the attachment point. The calling scope is
 * responsible for adding this direction to occupied_directions.
 * @throw NoAvailableDirectionsException if we run out of available directions.
 * This is guaranteed to not occur for at least the first eight attachment
 * points.
 */
static Direction calculate_direction_for_unbound_attachment_point(
    const int ap_num, const std::string& ap_name,
    const MonomerType monomer_type,
    const std::vector<BoundAttachmentPoint>& bound_aps,
    const std::vector<UnboundAttachmentPoint>& unbound_aps,
    const std::unordered_set<Direction>& occupied_directions)
{
    static const std::vector<Direction> DIAGONALS = {
        Direction::NW, Direction::NE, Direction::SE, Direction::SW};
    /**
     * @return the first available direction from the given list of cardinal
     * directions. If none of those directions are available, try the diagonals
     *
     * @throw NoAvailableDirectionsException if all of the given directions and
     * all of the diagonals are occupied
     */
    auto first_available = [&occupied_directions](
                               const std::vector<Direction>& dirs_to_try) {
        for (const auto cur_dir : boost::range::join(dirs_to_try, DIAGONALS)) {
            if (!occupied_directions.contains(cur_dir)) {
                return cur_dir;
            }
        }
        throw NoAvailableDirectionsException();
    };

    std::vector<Direction> dirs_to_try;
    if (ap_num <= 2 || ap_name == "pair") {
        // the first two attachment points should be across from each other (and
        // "pair" is the second attachment point for nucleic acid base monomers)
        int opposite_ap_num = ap_num == 2 || ap_name == "pair" ? 1 : 2;
        auto opposite_dir = fetch_assigned_direction_for_attachment_point(
            opposite_ap_num, bound_aps, unbound_aps);
        if (opposite_dir.has_value()) {
            // the other attachment point has a direction, so try to place this
            // attachment point on the opposite side
            auto preferred_dir = OPPOSITE_CARDINAL_DIRECTION.at(*opposite_dir);
            dirs_to_try.push_back(preferred_dir);
            // if there's already something on the opposite side, then we'll
            // have to put this somewhere else
            auto perpendiculars =
                PERPENDICULAR_CARDINAL_DIRECTIONS.at(*opposite_dir);
            std::copy(perpendiculars.begin(), perpendiculars.end(),
                      std::back_inserter(dirs_to_try));
            return first_available(dirs_to_try);
        } else if (monomer_type == MonomerType::NA_BASE) {
            // the other attachment point hasn't been laid out yet, which means
            // we must be laying out R1 in these next two conditions

            // prefer vertical for bases
            return first_available(
                {Direction::S, Direction::E, Direction::W, Direction::N});
        } else {
            // prefer horizontal for everything else
            return first_available(
                {Direction::W, Direction::N, Direction::S, Direction::E});
        }
    } else if (ap_num == 3) {
        // R3 (protein side chain interactions or the nucleic acid sugar to base
        // bond) should be perpendicular to R1 and R2
        auto r1_dir = fetch_assigned_direction_for_attachment_point(
            1, bound_aps, unbound_aps);
        if (!r1_dir.has_value()) {
            throw std::runtime_error("Directions for unbound attachment points "
                                     "must be calculated in order");
        }
        auto cw_dir = CLOCKWISE_CARDINAL_DIRECTION.at(*r1_dir);
        return first_available({cw_dir, OPPOSITE_CARDINAL_DIRECTION.at(cw_dir),
                                OPPOSITE_CARDINAL_DIRECTION.at(*r1_dir)});

    } else {
        // we should only get here with a CHEM monomer or a non-standard
        // attachment point (e.g. an R4 attachment point on an amino acid, which
        // normally only has three attachment points). Since we don't know what
        // this attachment point is supposed to represent, just find any
        // available direction.
        return first_available(
            {Direction::W, Direction::E, Direction::N, Direction::S});
    }
}

/**
 * @return A list of all attachment points on the specified monomer that are
 * currently unbound, with assigned directions for drawing the attachment point
 * "nubbin."
 *
 * Note that, if there are more attachment points than available directions, the
 * returned list will omit some attachment points if no available direction can
 * be assigned to them. This is guaranteed to only happen when there are more
 * than 8 attachment points for the monomer.
 *
 * Also note that we don't have a good way to determine how many attachment
 * points a CHEM monomer should have, so we assume that it has one additional
 * attachment point beyond the highest numbered bound attachment point.
 */
static std::vector<UnboundAttachmentPoint>
get_unbound_attachment_points(const RDKit::Atom* monomer,
                              std::vector<BoundAttachmentPoint> bound_aps)
{
    // build sets of bound attachment point numbers and custom names
    std::unordered_set<int> bound_ap_nums;
    std::unordered_set<std::string> bound_aps_with_custom_names;
    for (auto cur_ap : bound_aps) {
        if (cur_ap.num > 0) {
            bound_ap_nums.insert(cur_ap.num);
        } else if (cur_ap.num == ATTACHMENT_POINT_WITH_CUSTOM_NAME) {
            bound_aps_with_custom_names.insert(cur_ap.name);
        }
    }

    // figure out how many numbered attachment points we expect
    auto monomer_type = get_monomer_type(monomer);
    int num_numbered_aps = -1;
    if (NUMBERED_AP_NAMES_BY_MONOMER_TYPE.contains(monomer_type)) {
        num_numbered_aps =
            NUMBERED_AP_NAMES_BY_MONOMER_TYPE.at(monomer_type).size();
    } else if (monomer_type == MonomerType::NA_PHOSPHATE) {
        num_numbered_aps = 2;
    } else {
        // a CHEM monomer
        num_numbered_aps =
            *std::max_element(bound_ap_nums.begin(), bound_ap_nums.end());
        num_numbered_aps += 1;
    }

    std::unordered_set<Direction> occupied_directions;
    for (auto ap : bound_aps) {
        occupied_directions.insert(ap.direction);
    }

    std::vector<UnboundAttachmentPoint> available_aps;
    try {
        // figure out which numbered attachment points are unbound
        for (int ap_num = 1; ap_num <= num_numbered_aps; ++ap_num) {
            if (!bound_ap_nums.contains(ap_num)) {
                auto dir = calculate_direction_for_unbound_attachment_point(
                    ap_num, "", monomer_type, bound_aps, available_aps,
                    occupied_directions);
                occupied_directions.insert(dir);
                available_aps.push_back(
                    UnboundAttachmentPoint("", ap_num, dir));
            }
        }

        // figure out which attachment points with custom names are unbound
        if (EXPECTED_AP_CUSTOM_NAMES.contains(monomer_type)) {
            for (const auto& ap_name :
                 EXPECTED_AP_CUSTOM_NAMES.at(monomer_type)) {
                if (!bound_aps_with_custom_names.contains(ap_name)) {
                    auto dir = calculate_direction_for_unbound_attachment_point(
                        ATTACHMENT_POINT_WITH_CUSTOM_NAME, ap_name,
                        monomer_type, bound_aps, available_aps,
                        occupied_directions);
                    occupied_directions.insert(dir);
                    available_aps.push_back(UnboundAttachmentPoint(
                        ap_name, ATTACHMENT_POINT_WITH_CUSTOM_NAME, dir));
                }
            }
        }
    } catch (const NoAvailableDirectionsException&) {
        // there are so many attachment points that we can't generate directions
        // for all of them, so leave them out of the returned list
    }
    return available_aps;
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
            return ap_num_to_name(
                sugar_ap_num,
                NUMBERED_AP_NAMES_BY_MONOMER_TYPE.at(MonomerType::NA_SUGAR));
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

    if (NUMBERED_AP_NAMES_BY_MONOMER_TYPE.contains(monomer_type)) {
        return NUMBERED_AP_NAMES_BY_MONOMER_TYPE.at(monomer_type);
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

std::pair<std::vector<BoundAttachmentPoint>,
          std::vector<UnboundAttachmentPoint>>
get_attachment_points_for_monomer(const RDKit::Atom* monomer)
{
    auto bound_aps = get_bound_attachment_points(monomer);
    auto unbound_aps = get_unbound_attachment_points(monomer, bound_aps);

    // add pretty names to all numbered attachment points
    auto all_names = get_all_numbered_attachment_point_names(monomer);
    auto assign_ap_names = [&all_names](auto&& aps) {
        for (auto& cur_ap : aps) {
            if (cur_ap.num != ATTACHMENT_POINT_WITH_CUSTOM_NAME) {
                cur_ap.name = ap_num_to_name(cur_ap.num, all_names);
            }
        }
    };
    assign_ap_names(bound_aps);
    assign_ap_names(unbound_aps);

    return std::make_pair(bound_aps, unbound_aps);
}

std::string
get_attachment_point_name_for_connection(const RDKit::Atom* monomer,
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
