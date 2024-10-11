#include "schrodinger/sketcher/rdkit/rgroup.h"

#include <algorithm>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/RDGeneral/types.h>

#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/dummy_atom.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/molviewer/constants.h"

using schrodinger::rdkit_extensions::ATTACHMENT_POINT_LABEL_PREFIX;

namespace schrodinger
{
namespace sketcher
{

namespace
{
/**
 * Return the numerical suffix of a string.  E.g. if given arguments of ("_AP5",
 * "_AP"), return 5 as an integer.  If no numerical suffix is found, 0 will be
 * returned.
 *
 * @param label The string to fetch the suffix
 * @param prefix The prefix to strip from string
 * @return The numerical portion of the suffix
 */
unsigned int get_numerical_suffix(const std::string& label,
                                  const std::string& prefix)
{
    try {
        std::string digits = label.substr(prefix.size());
        return std::stoul(digits);
    } catch (const std::invalid_argument&) {
        return 0;
    } catch (const std::out_of_range&) {
        return 0;
    }
}
} // namespace

std::shared_ptr<RDKit::Atom>
make_new_attachment_point(const unsigned int ap_num)
{
    auto atom = rdkit_extensions::create_dummy_atom();
    atom->setProp(RDKit::common_properties::atomLabel,
                  ATTACHMENT_POINT_LABEL_PREFIX + std::to_string(ap_num));
    return atom;
}

bool is_attachment_point(const RDKit::Atom* const atom)
{
    // rdkit_extensions::is_attachment_point_dummy doesn't require digits after
    // "_AP", but we do since we need to use those numbers when we renumber the
    // attachment points after an atom deletion
    return get_attachment_point_number(atom);
}

unsigned int get_attachment_point_number(const RDKit::Atom* const atom)
{
    if (!rdkit_extensions::is_attachment_point_dummy(*atom)) {
        return 0;
    }
    // We know that the atomLabel property must be present, because
    // is_attachment_point_dummy will return false without it
    std::string label =
        atom->getProp<std::string>(RDKit::common_properties::atomLabel);
    return get_numerical_suffix(label, ATTACHMENT_POINT_LABEL_PREFIX);
}

bool is_r_group(const RDKit::Atom* const atom)
{
    return rdkit_extensions::get_r_group_number(atom).has_value();
}

std::vector<unsigned int> get_all_r_group_numbers(const RDKit::ROMol* const mol)
{
    std::unordered_set<unsigned int> r_groups;
    for (auto* atom : mol->atoms()) {
        if (auto r_group_num = rdkit_extensions::get_r_group_number(atom)) {
            r_groups.insert(r_group_num.value());
        }
    }
    std::vector<unsigned int> sorted_r_groups(r_groups.begin(), r_groups.end());
    std::sort(sorted_r_groups.begin(), sorted_r_groups.end());
    return sorted_r_groups;
}

std::vector<unsigned int>
get_next_r_group_numbers(const RDKit::ROMol* const mol, const size_t how_many)
{
    std::vector<unsigned int> sorted_r_groups = get_all_r_group_numbers(mol);
    std::vector<unsigned int> unused_r_groups;
    unused_r_groups.reserve(how_many);
    for (size_t i = 1; i <= sorted_r_groups.size(); ++i) {
        if (i != sorted_r_groups[i - 1]) {
            unused_r_groups.push_back(i);
            if (unused_r_groups.size() == how_many) {
                return unused_r_groups;
            }
            sorted_r_groups.insert(sorted_r_groups.begin() + i - 1, i);
        }
    }
    unsigned int next_r =
        sorted_r_groups.empty() ? 1 : sorted_r_groups.back() + 1;
    for (unsigned int i = next_r; unused_r_groups.size() < how_many; ++i) {
        unused_r_groups.push_back(i);
    }
    return unused_r_groups;
}

unsigned int get_next_attachment_point_number(const RDKit::ROMol* const mol)
{
    unsigned int max_ap_num = 0;
    for (auto* atom : mol->atoms()) {
        max_ap_num = std::max(max_ap_num, get_attachment_point_number(atom));
    }
    return max_ap_num + 1;
}

void renumber_attachment_points(RDKit::RWMol* const mol)
{
    std::vector<RDKit::Atom*> attachment_points;
    for (auto* atom : mol->atoms()) {
        if (is_attachment_point(atom)) {
            attachment_points.push_back(atom);
        }
    }
    std::sort(attachment_points.begin(), attachment_points.end(),
              [](RDKit::Atom* a, RDKit::Atom* b) {
                  return get_attachment_point_number(a) <
                         get_attachment_point_number(b);
              });
    unsigned int ap_num = 1;
    for (auto* atom : attachment_points) {
        atom->setProp(RDKit::common_properties::atomLabel,
                      ATTACHMENT_POINT_LABEL_PREFIX + std::to_string(ap_num++));
    }
}

const RDKit::Bond* const
get_attachment_point_bond(const RDKit::Atom* const atom)
{
    if (is_attachment_point(atom)) {
        auto& mol = atom->getOwningMol();
        return *(mol.atomBonds(atom).begin());
    }
    return nullptr;
}

const RDKit::Atom* const
get_attachment_point_atom(const RDKit::Bond* const bond)
{
    const auto* begin_atom = bond->getBeginAtom();
    const auto* end_atom = bond->getEndAtom();
    if (is_attachment_point(begin_atom)) {
        return begin_atom;
    } else if (is_attachment_point(end_atom)) {
        return end_atom;
    }
    return nullptr;
}

bool is_attachment_point_bond(const RDKit::Bond* const bond)
{
    return get_attachment_point_atom(bond) != nullptr;
}

unsigned int number_of_bound_attachment_points(const RDKit::Atom* const atom)
{
    auto& mol = atom->getOwningMol();
    auto neighbors = mol.atomNeighbors(atom);
    return std::count_if(neighbors.begin(), neighbors.end(),
                         [](const RDKit::Atom* const cur_neighbor) {
                             return is_attachment_point(cur_neighbor);
                         });
}

} // namespace sketcher
} // namespace schrodinger
