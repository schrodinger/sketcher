#include "schrodinger/sketcher/rdkit/variable_attachment_bond_core.h"

#include <string>
#include <sstream>

#include <boost/algorithm/string/trim.hpp>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/Depictor/DepictUtils.h>

#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/rdkit_extensions/dummy_atom.h"
#include "schrodinger/sketcher/rdkit/coord_utils.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * The value for the _MolFileBondAttach property (i.e. the ATTACH value in an
 * MDL molfile) of newly-created variable attachment bonds
 */
const std::string MOL_FILE_BOND_ATTACH = "ANY";

bool is_variable_attachment_bond(const RDKit::Bond* bond)
{
    return !get_variable_attachment_atoms(bond).empty();
}

bool is_dummy_atom_for_variable_attachment_bond(const RDKit::Atom* atom)
{
    // is this a dummy atom with exactly one bond that's a variable attachment
    // bond?
    return (atom->getAtomicNum() == rdkit_extensions::DUMMY_ATOMIC_NUMBER &&
            atom->getDegree() == 1 &&
            is_variable_attachment_bond(
                *(atom->getOwningMol().atomBonds(atom).begin())));
}

std::unordered_set<const RDKit::Atom*>
get_variable_attachment_atoms(const RDKit::Bond* bond)
{
    std::string end_pts_text;
    bond->getPropIfPresent(RDKit::common_properties::_MolFileBondEndPts,
                           end_pts_text);
    if (end_pts_text.empty()) {
        // the bond is missing the end points property, or it's empty
        return {};
    }
    bool begin_atom_is_dummy = bond->getBeginAtom()->getAtomicNum() ==
                               rdkit_extensions::DUMMY_ATOMIC_NUMBER;
    bool end_atom_is_dummy = bond->getEndAtom()->getAtomicNum() ==
                             rdkit_extensions::DUMMY_ATOMIC_NUMBER;
    if (begin_atom_is_dummy == end_atom_is_dummy) {
        // the bond must have exactly one bound dummy atom
        return {};
    }
    boost::trim(end_pts_text);
    if (!(end_pts_text.front() == '(' && end_pts_text.back() == ')')) {
        // the end points property must start and end with parenthesis
        return {};
    }
    // remove the parenthesis and any leading or trailing whitespace
    end_pts_text.pop_back();
    end_pts_text.erase(end_pts_text.begin());
    boost::trim(end_pts_text);

    // convert the string to a set of valid atom indices
    const auto& mol = bond->getOwningMol();
    auto num_atoms = static_cast<int>(mol.getNumAtoms());
    std::unordered_set<unsigned int> atom_indices;
    std::istringstream ss(end_pts_text);
    int num_bound_atoms;
    ss >> num_bound_atoms;
    if (ss.fail() || num_bound_atoms < 0) {
        return {};
    }
    while (!ss.eof()) {
        int cur_index;
        ss >> cur_index;
        if (ss.fail() || cur_index < 1 || cur_index > num_atoms) {
            return {};
        }
        // MDL molfiles are 1-indexed instead of 0-indexed, so we have to
        // subtract one here
        atom_indices.insert(static_cast<unsigned int>(cur_index) - 1);
    }

    // convert the atom indices to atoms
    std::unordered_set<const RDKit::Atom*> atoms;
    std::transform(atom_indices.begin(), atom_indices.end(),
                   std::inserter(atoms, atoms.begin()),
                   [&mol](const auto idx) { return mol.getAtomWithIdx(idx); });
    return atoms;
}

/**
 * Generate the set of potential crossing bonds for a variable attachment bond.
 * The crossing bond controls where we draw the variable attachment bond, which
 * is always perpendicular to the crossing bond.
 * @param mol The molecule we are adding the variable attachment bond to
 * @param atoms The variable attachment atoms
 * @return The potential crossing bonds. These are all bonds between the
 * variable attachment atoms. If no such bonds exist, then the potential
 * crossing bonds are any bonds involving the variable attachment atoms
 * @throw std::runtime_error if none of the variable attachment atoms are
 * involved in any bonds
 */
static std::unordered_set<const RDKit::Bond*> get_potential_crossing_bonds(
    const RDKit::ROMol& mol,
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    std::unordered_set<const RDKit::Bond*> bonds_between_specified_atoms;
    std::unordered_set<const RDKit::Bond*> bonds_involving_specified_atoms;
    for (auto* cur_atom : atoms) {
        for (auto* cur_bond : mol.atomBonds(cur_atom)) {
            bonds_involving_specified_atoms.insert(cur_bond);
            if (atoms.count(cur_bond->getOtherAtom(cur_atom))) {
                bonds_between_specified_atoms.insert(cur_bond);
            }
        }
    }
    if (!bonds_between_specified_atoms.empty()) {
        return bonds_between_specified_atoms;
    } else if (!bonds_involving_specified_atoms.empty()) {
        return bonds_involving_specified_atoms;
    } else {
        // we're guaranteed to never hit this if all variable attachment atoms
        // are in the same molecule, since that ensures that the atoms have at
        // least one bond
        throw variable_attachment_bond_error("No bonds found");
    }
}

/**
 * @return The midpoint of the specified bond
 */
static RDGeom::Point3D get_bond_midpoint(const RDKit::Bond* bond)
{
    const auto& conf = bond->getOwningMol().getConformer();
    const auto& start_coords = conf.getAtomPos(bond->getBeginAtomIdx());
    const auto& end_coords = conf.getAtomPos(bond->getEndAtomIdx());
    return (start_coords + end_coords) / 2.0;
}

/**
 * Pick the best crossing bond from a list of potential crossing bonds. The
 * crossing bond controls where we draw the variable attachment bond, which is
 * always perpendicular to the crossing bond. The best crossing bond is the one
 * with the minimum distance from the crossing bond's midpoint to the centroid
 * of all variable attachment atoms.
 * @param potential_crossing_bonds All potential crossing bonds, as generated by
 * get_potential_crossing_bonds
 * @param centroid The centroid of all variable attachment atoms
 * @return The crossing bond
 */
static const RDKit::Bond* get_best_crossing_bond(
    const std::unordered_set<const RDKit::Bond*>& potential_crossing_bonds,
    const RDGeom::Point3D& centroid)
{
    auto dist_sq_from_bond_midpoint_to_centroid =
        [&centroid](const RDKit::Bond* bond) {
            auto midpoint = get_bond_midpoint(bond);
            return (centroid - midpoint).lengthSq();
        };
    return *std::min_element(
        potential_crossing_bonds.begin(), potential_crossing_bonds.end(),
        [&dist_sq_from_bond_midpoint_to_centroid](auto* bond_a, auto* bond_b) {
            return dist_sq_from_bond_midpoint_to_centroid(bond_a) <
                   dist_sq_from_bond_midpoint_to_centroid(bond_b);
        });
}

/**
 * Generate the coordinates for the ends of a variable attachment bond. The
 * variable attachment bond will be 1.5x as long as a regular bond, with 1/3 of
 * that length inside the ring (i.e. the dummy atom side of the crossing bond)
 * and 2/3 outside the ring (i.e. the carbon atom side of the crossing bond).
 * @param mol The molecule we are adding the variable attachment bond to
 * @param atoms The variable attachment atoms
 * @param crossing_bond The bond to draw the variable attachment bond
 * perpendicular to
 * @return A pair of
 *   - the coordinates for the variable attachment bond's dummy atom
 *   - the coordinates for the variable attachment bond's carbon atom
 */
static std::pair<RDGeom::Point3D, RDGeom::Point3D>
get_variable_bond_endpoints(const RDKit::ROMol& mol,
                            const std::unordered_set<const RDKit::Atom*>& atoms,
                            const RDKit::Bond* crossing_bond)
{
    auto midpoint = get_bond_midpoint(crossing_bond);

    const auto& conf = crossing_bond->getOwningMol().getConformer();
    const auto& start_coords =
        conf.getAtomPos(crossing_bond->getBeginAtomIdx());
    auto normal = (midpoint - start_coords).getPerpendicular();
    normal *= (RDDepict::BOND_LEN / 2.0);
    auto dummy_coords_a = midpoint + normal;
    auto dummy_coords_b = midpoint - normal;

    // Get all the neighbors of the variable attachment atoms that aren't
    // variable attachment atoms themselves
    std::unordered_set<const RDKit::Atom*> neighbors;
    for (const auto* cur_atom : atoms) {
        for (const auto* cur_neighbor : mol.atomNeighbors(cur_atom)) {
            if (!atoms.count(cur_neighbor)) {
                neighbors.insert(cur_neighbor);
            }
        }
    }
    auto closest_neighbor_dist_sq =
        [&mol, &neighbors](const RDGeom::Point3D& dummy_coords) {
            const auto& conf = mol.getConformer();
            std::vector<double> dist_sqs;
            std::transform(neighbors.begin(), neighbors.end(),
                           std::back_inserter(dist_sqs),
                           [&conf, &dummy_coords](auto* atom) {
                               const auto& atom_pos =
                                   conf.getAtomPos(atom->getIdx());
                               return (atom_pos - dummy_coords).lengthSq();
                           });
            return *std::min_element(dist_sqs.begin(), dist_sqs.end());
        };

    RDGeom::Point3D dummy_coords;
    RDGeom::Point3D real_atom_coords;
    normal *= 2;
    if (closest_neighbor_dist_sq(dummy_coords_a) <=
        closest_neighbor_dist_sq(dummy_coords_b)) {
        dummy_coords = dummy_coords_a;
        real_atom_coords = midpoint - normal;
    } else {
        dummy_coords = dummy_coords_b;
        real_atom_coords = midpoint + normal;
    }
    return {dummy_coords, real_atom_coords};
}

/**
 * Generate the text for the variable attachment bond's ENDPTS property
 * @param atoms The variable attachment atoms
 * @return The ENDPTS text
 */
static std::string generate_bond_end_pts_property_text(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    std::vector<unsigned int> atom_idxs;
    // MDL molfiles are 1-indexed instead of 0-indexed, so we have to add one
    // here
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(atom_idxs),
                   [](auto* atom) { return atom->getIdx() + 1; });
    std::sort(atom_idxs.begin(), atom_idxs.end());
    std::ostringstream ss;
    ss << "(" << atom_idxs.size() << " ";
    std::copy(atom_idxs.begin(), atom_idxs.end(),
              std::ostream_iterator<unsigned int>(ss, " "));
    // replace the last space with a close parenthesis
    ss.seekp(-1, std::ios_base::end);
    ss << ")";
    return ss.str();
}

/**
 * Create the dummy atom, carbon atom, and variable attachment bond, and then
 * add them to the molecule
 * @param mol The molecule to add the atoms and bond to
 * @param dummy_coords The coordinates for the dummy atom
 * @param real_atom_coords The coordinates for the carbon atom
 * @param bond_end_pts_text The text of the variable attachment bond's ENDPTS
 * property
 * @return A tuple of
 *   - the dummy atom
 *   - the carbon atom
 *   - the variable attachment bond
 */
static std::tuple<RDKit::Atom*, RDKit::Atom*, RDKit::Bond*>
create_and_add_atoms_and_variable_attachment_bond(
    RDKit::RWMol& mol, const RDGeom::Point3D& dummy_coords,
    const RDGeom::Point3D& real_atom_coords, std::string bond_end_pts_text)
{
    auto dummy_atom = rdkit_extensions::create_dummy_atom();
    auto carbon = std::make_shared<RDKit::Atom>("C");
    auto bond = std::make_shared<RDKit::Bond>(RDKit::Bond::BondType::SINGLE);
    bond->setProp(RDKit::common_properties::_MolFileBondAttach,
                  MOL_FILE_BOND_ATTACH);
    bond->setProp(RDKit::common_properties::_MolFileBondEndPts,
                  bond_end_pts_text);

    unsigned int dummy_idx =
        mol.addAtom(dummy_atom.get(), /* updateLabel = */ false);
    unsigned int carbon_idx =
        mol.addAtom(carbon.get(), /* updateLabel = */ false);
    auto& conf = mol.getConformer();
    conf.setAtomPos(dummy_idx, dummy_coords);
    conf.setAtomPos(carbon_idx, real_atom_coords);

    bond->setOwningMol(mol);
    bond->setBeginAtomIdx(dummy_idx);
    bond->setEndAtomIdx(carbon_idx);
    unsigned int bond_idx = mol.addBond(bond.get()) - 1;

    return {mol.getAtomWithIdx(dummy_idx), mol.getAtomWithIdx(carbon_idx),
            mol.getBondWithIdx(bond_idx)};
}

std::pair<RDGeom::Point3D, RDGeom::Point3D>
get_coordinates_for_variable_attachment_bond(
    const RDKit::ROMol& mol,
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    if (atoms.size() < 2) {
        throw variable_attachment_bond_error(
            "Variable attachment bonds require at least two attachment atoms.");
    }
    auto centroid = find_centroid(mol, atoms);
    auto potential_crossing_bonds = get_potential_crossing_bonds(mol, atoms);
    auto crossing_bond =
        get_best_crossing_bond(potential_crossing_bonds, centroid);
    return get_variable_bond_endpoints(mol, atoms, crossing_bond);
}

std::tuple<RDKit::Atom*, RDKit::Atom*, RDKit::Bond*>
add_variable_attachment_bond_to_mol(
    RDKit::RWMol& mol, const std::unordered_set<const RDKit::Atom*>& atoms)
{
    auto [dummy_coords, real_atom_coords] =
        get_coordinates_for_variable_attachment_bond(mol, atoms);
    auto bond_end_pts_text = generate_bond_end_pts_property_text(atoms);
    return create_and_add_atoms_and_variable_attachment_bond(
        mol, dummy_coords, real_atom_coords, bond_end_pts_text);
}

} // namespace sketcher
} // namespace schrodinger
