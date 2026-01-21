#include "schrodinger/sketcher/rdkit/monomeric.h"

#include <functional>

#include <fmt/core.h>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>

#include "schrodinger/rdkit_extensions/helm.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{
const std::string PEPTIDE_POLYMER_PREFIX = "PEPTIDE";
// According to HELM, DNA is a subtype of RNA, so DNA also uses the RNA prefix
const std::string NUCLEOTIDE_POLYMER_PREFIX = "RNA";

// note that the primes use apostrophes instead of a Unicode prime to avoid
// issues with C++′s handling of Unicode
const std::unordered_map<MonomerType, std::vector<std::string>> AP_NAMES = {
    {MonomerType::PEPTIDE, {"N", "C", "X"}},
    {MonomerType::NA_BASE, {"S", "BP"}},
    {MonomerType::NA_SUGAR, {"3'", "5'", "X"}},
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

bool contains_two_monomer_linkages(const RDKit::Bond* bond)
{
    std::string linkage, custom_linkage;
    bond->getPropIfPresent(LINKAGE, linkage);
    bool custom_linkage_exists =
        bond->getPropIfPresent(CUSTOM_BOND, custom_linkage);
    return custom_linkage_exists && custom_linkage != linkage;
}

/**
 * Return the number of the attachment point specified in the given linkage
 * string.
 * @param linkage A description of the linkage taken from the bond properties.
 * Formatted similar to "R2-R3".
 * @param is_begin_atom If true, we'll return the attachment point for the
 * bond's begin atom.  Otherwise, we'll return the attachment point for the
 * bond's end atom.
 * @return the attachment point number, or -1 if the attachment point is not
 * properly specified
 */
static int get_attachment_point_for_atom(std::string linkage,
                                         bool is_begin_atom)
{
    auto dash_pos = linkage.find("-");
    if (dash_pos == std::string::npos) {
        return -1;
    }
    std::string attachment_point_name = is_begin_atom
                                            ? linkage.substr(0, dash_pos)
                                            : linkage.substr(dash_pos + 1);
    if (attachment_point_name[0] != 'R') {
        return -1;
    }
    // remove the leading 'R' now that we've confirmed it exists
    attachment_point_name.erase(0, 1);
    try {
        return std::stoi(attachment_point_name);
    } catch (const std::logic_error&) {
        // it's not an integer
        return -1;
    }
}

/**
 * @overload Return the number of the attachment point for the bond on the
 * specified monomer. Note that if the bond specifies two linkages, this
 * overload will only return the attachment point of the primary linkage.
 * @return the attachment point number, or -1 if the attachment point is not
 * properly specified
 */
static int get_attachment_point_for_atom(const RDKit::Atom* monomer,
                                         const RDKit::Bond* bond)
{
    std::string linkage;
    if (bond->getPropIfPresent(LINKAGE, linkage)) {
        bool is_start_atom = bond->getBeginAtom() == monomer;
        return get_attachment_point_for_atom(linkage, is_start_atom);
    }
    return -1;
}

/**
 * @return a map of {attachment point: bound monomer} for all bound attachment
 * points of the specified monomer. Attachment points are specified using
 * integers, e.g. 1 for "R1".
 */
static std::vector<std::pair<int, const RDKit::Atom*>>
get_bound_attachment_points(const RDKit::Atom* monomer)
{
    const auto& mol = monomer->getOwningMol();
    std::unordered_set<int> bound_ap_nums;
    std::vector<std::pair<int, const RDKit::Atom*>> bound_aps;

    auto record_linkage = [&bound_aps, &bound_ap_nums, monomer](
                              const RDKit::Bond* bond, const bool is_start_atom,
                              const std::string& prop_name) {
        std::string linkage;
        if (bond->getPropIfPresent(prop_name, linkage)) {
            auto ap_value =
                get_attachment_point_for_atom(linkage, is_start_atom);
            if (ap_value > 0 && !bound_ap_nums.contains(ap_value)) {
                bound_ap_nums.insert(ap_value);
                bound_aps.push_back({ap_value, bond->getOtherAtom(monomer)});
            }
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
 * @return a set of all attachment points on the specified monomer that are
 * currently unbound. Attachment points are specified using integers, e.g. 1 for
 * "R1".
 *
 * Note that we don't have a good way to determine how many attachment points a
 * CHEM monomer should have, so we assume that it has one additional attachment
 * point beyond the highest numbered bound attachment point.
 */
static std::vector<int>
get_available_attachment_points(const RDKit::Atom* monomer)
{
    auto bound_aps = get_bound_attachment_points(monomer);
    auto monomer_type = get_monomer_type(monomer);
    std::unordered_set<int> bound_ap_nums;
    std::transform(bound_aps.begin(), bound_aps.end(),
                   std::inserter(bound_ap_nums, bound_ap_nums.end()),
                   [](auto num_and_atom) { return num_and_atom.first; });
    int num_aps = -1;
    if (AP_NAMES.contains(monomer_type)) {
        num_aps = AP_NAMES.at(monomer_type).size();
    } else if (monomer_type == MonomerType::NA_PHOSPHATE) {
        num_aps = 2;
    } else {
        // a CHEM monomer
        num_aps = *std::max_element(bound_ap_nums.begin(), bound_ap_nums.end());
        num_aps += 1;
    }
    std::vector<int> available_aps;
    for (int ap = 1; ap <= num_aps; ++ap) {
        if (!bound_ap_nums.contains(ap)) {
            available_aps.push_back(ap);
        }
    }
    return available_aps;
}

/**
 * Convert an attachment point number to a name
 * @param ap_num The attachment point number to convert
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
        auto sugar_ap_num =
            get_attachment_point_for_atom(cur_neighbor, bond_to_sugar);
        // the phosphate should be bound to either the 3' (R1) or 5' (R2). If
        // it's bound to something else, ignore it since something's gone wrong.
        if ((sugar_ap_num == 1 || sugar_ap_num == 2)) {
            return ap_num_to_name(sugar_ap_num,
                                  AP_NAMES.at(MonomerType::NA_SUGAR));
        }
    }
    return "";
}

/**
 * @return a list of all "pretty" attachment point names (e.g. "N" instead of
 * "R1" for amino acids) for the given monomer, regardless of whether those
 * attachment points are bound or available.
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
get_all_attachment_point_names(const RDKit::Atom* monomer)
{
    std::vector<std::string> all_names;
    auto monomer_type = get_monomer_type(monomer);

    if (AP_NAMES.contains(monomer_type)) {
        return AP_NAMES.at(monomer_type);
    } else if (monomer_type == MonomerType::NA_PHOSPHATE) {
        std::vector<std::string> phos_ap_names = {"", ""};
        auto sugar_ap_name = get_attachment_point_name_of_bound_sugar(monomer);
        if (!sugar_ap_name.empty()) {
            const auto& mol = monomer->getOwningMol();
            // the phosphate must have exactly one bond; otherwise,
            // sugar_ap_name would be empty
            const RDKit::Bond* phos_bond = *mol.atomBonds(monomer).begin();
            auto bound_phos_ap_num =
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
    return all_names;
}

std::vector<std::pair<std::string, const RDKit::Atom*>>
get_bound_attachment_point_names_and_atoms(const RDKit::Atom* monomer)
{
    auto bound_aps = get_bound_attachment_points(monomer);
    auto all_names = get_all_attachment_point_names(monomer);
    std::vector<std::pair<std::string, const RDKit::Atom*>> bound_ap_names;
    std::transform(bound_aps.begin(), bound_aps.end(),
                   std::back_inserter(bound_ap_names),
                   [&all_names](auto num_and_atom) {
                       auto& [ap_num, atom] = num_and_atom;
                       auto ap_name = ap_num_to_name(ap_num, all_names);
                       return std::make_pair(ap_name, atom);
                   });
    return bound_ap_names;
}

std::vector<std::string>
get_available_attachment_point_names(const RDKit::Atom* monomer)
{
    auto available_aps = get_available_attachment_points(monomer);
    auto all_names = get_all_attachment_point_names(monomer);
    std::vector<std::string> available_names;
    std::transform(available_aps.begin(), available_aps.end(),
                   std::back_inserter(available_names),
                   std::bind(ap_num_to_name, std::placeholders::_1, all_names));
    return available_names;
}

} // namespace sketcher
} // namespace schrodinger
