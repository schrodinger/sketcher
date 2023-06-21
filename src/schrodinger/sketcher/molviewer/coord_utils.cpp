#include "schrodinger/sketcher/molviewer/coord_utils.h"

#include <algorithm>
#include <iterator>

#include <boost/range/iterator_range.hpp>

#include <Geometry/point.h>
#include <GraphMol/CoordGen.h>
#include <GraphMol/ROMol.h>

#include <QLineF>

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

SKETCHER_API bool has_2d_conformer(const RDKit::ROMol& mol)
{
    auto conformer =
        std::find_if_not(mol.beginConformers(), mol.endConformers(),
                         std::mem_fn(&RDKit::Conformer::is3D));

    return conformer != mol.endConformers();
}

SKETCHER_API RDKit::Conformer& get_2d_conformer(RDKit::ROMol& mol)
{
    auto conformer =
        std::find_if_not(mol.beginConformers(), mol.endConformers(),
                         std::mem_fn(&RDKit::Conformer::is3D));

    if (conformer != mol.endConformers()) {
        return *(*conformer);
    }

    const auto conformer_id = RDKit::CoordGen::addCoords(mol);
    return mol.getConformer(conformer_id);
}

SKETCHER_API void rescale(RDKit::Conformer& conformer, RDKit::ROMol& mol,
                          double reference_bond_length)
{
    if (&conformer.getOwningMol() != &mol) {
        throw std::invalid_argument(
            "Conformer does not belong to the given molecule");
    }

    if (mol.getNumBonds() == 0) {
        return;
    }

    const double tolerance = 0.1;

    // build a map storing how many times each bond length is present
    std::unordered_map<double, int> length_map;
    for (auto bond : mol.bonds()) {
        auto this_bond_length = (conformer.getAtomPos(bond->getBeginAtomIdx()),
                                 conformer.getAtomPos(bond->getEndAtomIdx()))
                                    .length();

        // check if we have already seen a similar bond length
        for (auto& length_pair : length_map) {
            double ratio = length_pair.first / this_bond_length;
            if (ratio < 1. + tolerance && ratio > 1. - tolerance) {
                this_bond_length = length_pair.first;
            }
            break;
        }

        // Make sure we don't end up dividing by 0 (or something tiny)
        if (std::abs(this_bond_length) > 0.0001) {
            length_map[this_bond_length]++;
        }
    }

    // Make sure we at least have one bond to measure
    if (length_map.empty()) {
        return;
    }

    auto most_frequent_pair =
        std::max_element(length_map.begin(), length_map.end(),
                         [](const auto left, const auto& right) {
                             return left.second < right.second;
                         });
    auto ratio = reference_bond_length / most_frequent_pair->first;
    for (auto& atom_position : conformer.getPositions()) {
        atom_position *= ratio;
    }
}

QPointF to_scene_xy(const RDGeom::Point3D& mol_xy)
{
    return QPointF(mol_xy.x * VIEW_SCALE, -mol_xy.y * VIEW_SCALE);
}

RDGeom::Point3D to_mol_xy(const QPointF& scene_xy)
{
    return RDGeom::Point3D(scene_xy.x() / VIEW_SCALE,
                           -scene_xy.y() / VIEW_SCALE, 0.0);
}

std::vector<RDGeom::Point3D>
get_relative_positions_of_atom_neighbors(const RDKit::Atom* const atom)
{
    auto& mol = atom->getOwningMol();
    auto& conf = mol.getConformer();
    std::vector<RDGeom::Point3D> positions;
    const RDGeom::Point3D& atom_pos = conf.getAtomPos(atom->getIdx());

    for (const auto& neighbor : mol.atomNeighbors(atom)) {
        const RDGeom::Point3D& neighbor_pos =
            conf.getAtomPos(neighbor->getIdx());
        positions.push_back(neighbor_pos - atom_pos);
    }

    return positions;
}

QPointF best_placing_around_origin(const std::vector<QPointF>& points)
{
    if (points.empty()) {
        return QPointF(VIEW_SCALE, 0);
    }
    QPointF origin(0.f, 0.f);
    std::vector<float> angles;
    // find out the angles at which each member of points is found around the
    // origin
    for (auto& point : points) {
        QLineF line(origin, point);
        auto angle = line.angle();
        angles.push_back(angle);
    }
    sort(angles.begin(), angles.end());
    angles.push_back(angles.front() + 360);
    // find the biggest angle interval and return a point in the middle of it
    int best_i = 0;
    auto best_angle = (angles[best_i + 1] - angles[best_i]) * 0.5;
    for (unsigned int i = 0; i < angles.size() - 1; ++i) {
        auto angle_i = (angles[i + 1] - angles[i]) * 0.5;
        if (angle_i > best_angle) {
            best_i = i;
            best_angle = angle_i;
        }
    }
    // For a single substituent, limit the angle to 120 (instead of 180)
    float max_angle = (points.size() == 1) ? 120.0 : 180.0;
    float angle_to_use = (angles[best_i + 1] - angles[best_i]) * 0.5;
    float angle_increment = std::min(max_angle, angle_to_use);
    float return_angle = angles[best_i] + angle_increment;
    QLineF line(origin, QPointF(VIEW_SCALE, 0));
    line.setAngle(return_angle);
    return line.p2();
}

RDGeom::Point3D
best_placing_around_origin(const std::vector<RDGeom::Point3D>& points)
{
    std::vector<QPointF> qpoints;
    qpoints.reserve(points.size());
    std::transform(points.begin(), points.end(), std::back_inserter(qpoints),
                   to_scene_xy);
    QPointF best = best_placing_around_origin(qpoints);
    return to_mol_xy(best);
}

} // namespace sketcher
} // namespace schrodinger
