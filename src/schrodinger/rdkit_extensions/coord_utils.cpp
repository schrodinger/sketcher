#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"

#include <algorithm>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/Depictor/DepictUtils.h>
#include <rdkit/GraphMol/Depictor/RDDepictor.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{

/**
 * The minimum distance between two points before this module considers them as
 * separate coordinates.
 */
constexpr double BOND_LENGTH_EPSILON = 0.001;
constexpr double BOND_LENGTH_EPSILON_SQ =
    BOND_LENGTH_EPSILON * BOND_LENGTH_EPSILON;

/**
 * Check if the coordinates in the given conformer are all zero.  This can
 * happen when exporting an RDKit molecule without coordinates to a format that
 * requires coordinates, e.g., PDB.
 */
bool coordinates_are_all_zero(const RDKit::Conformer& conf)
{
    auto positions = conf.getPositions();
    return std::all_of(positions.begin(), positions.end(), [](auto pos) {
        return pos.lengthSq() < BOND_LENGTH_EPSILON_SQ;
    });
}

} // namespace

unsigned int compute2DCoords(RDKit::ROMol& mol,
                             const std::vector<unsigned int>& frozen_ids)
{
    if (isMonomeric(mol)) {
        return compute_monomer_mol_coords(mol);
    }

    RDDepict::Compute2DCoordParameters params;
    params.forceRDKit = true;
    params.useRingTemplates = true;

    RDGeom::INT_POINT2D_MAP coord_map;
    if (!frozen_ids.empty()) {
        auto& conf = mol.getConformer();
        for (auto idx : frozen_ids) {
            coord_map[idx] = conf.getAtomPos(idx);
        }
        params.coordMap = &coord_map;
    }

    // If we are generating a new 2D conformer, we don't want
    // to keep the old bond wedging, as it is relative to the old
    // conformation, which is overwritten by the new one.
    for (auto bond : mol.bonds()) {
        if (bond->hasProp(RDKit::common_properties::_MolFileBondStereo)) {
            bond->clearProp(RDKit::common_properties::_MolFileBondStereo);
            continue;
        }
        if (bond->hasProp(RDKit::common_properties::_MolFileBondCfg)) {
            bond->clearProp(RDKit::common_properties::_MolFileBondCfg);
        }
    }

    return RDDepict::compute2DCoords(mol, params);
}

void update_2d_coordinates(RDKit::ROMol& mol)
{

    ::RDKit::Conformer* conf_to_keep_3d = nullptr;
    auto conf_3d = std::find_if(mol.beginConformers(), mol.endConformers(),
                                std::mem_fn(&RDKit::Conformer::is3D));
    if (conf_3d != mol.endConformers() &&
        !coordinates_are_all_zero(**conf_3d)) {
        conf_to_keep_3d = new ::RDKit::Conformer(**conf_3d);
    }

    auto conf_2d = std::find_if_not(mol.beginConformers(), mol.endConformers(),
                                    std::mem_fn(&RDKit::Conformer::is3D));
    if (conf_2d == mol.endConformers() || coordinates_are_all_zero(**conf_2d)) {
        // remove all (3d) conformers and generate new 2d coordinates
        mol.clearConformers();
        compute2DCoords(mol);
    } else {
        // keep only the first 2D conformer and rescale the bond length
        auto conf_to_keep_2d = new ::RDKit::Conformer(**conf_2d);
        mol.clearConformers();
        mol.addConformer(conf_to_keep_2d, true);
        rescale_bond_length_if_needed(mol);
    }

    if (conf_to_keep_3d != nullptr) {
        mol.addConformer(conf_to_keep_3d, true);
    }
}

double get_most_common_bond_length(const RDKit::ROMol& mol)
{
    if (!mol.getNumConformers() || !mol.getNumBonds()) {
        return -1.0;
    }
    const auto& conf = mol.getConformer();
    std::unordered_map<int, std::vector<double>> bond_lengths_by_nearest_tenth;
    for (auto* cur_bond : mol.bonds()) {
        auto begin_pos = conf.getAtomPos(cur_bond->getBeginAtomIdx());
        auto end_pos = conf.getAtomPos(cur_bond->getEndAtomIdx());
        auto length = (end_pos - begin_pos).length();
        int nearest_tenth = std::round(length * 10);
        bond_lengths_by_nearest_tenth.try_emplace(nearest_tenth);
        bond_lengths_by_nearest_tenth[nearest_tenth].push_back(length);
    }
    auto most_common_length_it = std::max_element(
        bond_lengths_by_nearest_tenth.begin(),
        bond_lengths_by_nearest_tenth.end(),
        [](auto a, auto b) { return a.second.size() < b.second.size(); });
    auto lengths = most_common_length_it->second;
    auto average_length =
        std::reduce(lengths.begin(), lengths.end()) / lengths.size();
    return average_length;
}

void rescale_bond_length_if_needed(RDKit::ROMol& mol)
{
    auto cur_bond_length = get_most_common_bond_length(mol);
    if (cur_bond_length < BOND_LENGTH_EPSILON) {
        // If cur_bond_length is -1, then there's no conformer or no bonds in
        // the molecule, so there's nothing to scale.  If the bond length is
        // zero or near-zero, then the molecule is really messed up and we're
        // going to run into underflow or overflow issues, so give up.
        return;
    } else if (std::abs(RDDepict::BOND_LEN - cur_bond_length) <
               BOND_LENGTH_EPSILON) {
        // The bond length is the same as the expected bond length, so we don't
        // need to do anything
        return;
    }
    double scale = RDDepict::BOND_LEN / cur_bond_length;
    for (auto& coord : mol.getConformer().getPositions()) {
        coord *= scale;
    }
}

RDGeom::Point3D
find_centroid(const RDKit::ROMol& mol,
              const std::unordered_set<const RDKit::Atom*>& atoms)

{
    return find_centroid(mol.getConformer(), atoms);
}

RDGeom::Point3D
find_centroid(const RDKit::Conformer& conf,
              const std::unordered_set<const RDKit::Atom*>& atoms)
{
    std::vector<RDGeom::Point3D> positions;
    if (atoms.empty()) {
        positions = conf.getPositions();
    } else {
        positions.reserve(atoms.size());
        for (const auto& atom : atoms) {
            positions.push_back(conf.getAtomPos(atom->getIdx()));
        }
    }
    return find_centroid(positions);
}

RDGeom::Point3D find_centroid(const std::vector<RDGeom::Point3D>& positions)
{
    // Calculate the centroid by averaging the coordinates
    size_t num_atoms = positions.size();
    RDGeom::Point3D centroid;
    for (const auto& coord : positions) {
        centroid += coord;
    }
    // avoid division by zero
    if (num_atoms > 0) {
        centroid /= static_cast<double>(num_atoms);
    }
    return centroid;
}

} // namespace rdkit_extensions
} // namespace schrodinger
