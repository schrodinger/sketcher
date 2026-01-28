#include "schrodinger/sketcher/molviewer/coord_utils.h"

#include <algorithm>
#include <cmath>
#include <cfloat>
#include <iterator>

#include <boost/range/iterator_range.hpp>

#include <rdkit/Geometry/point.h>
#include <rdkit/GraphMol/ROMol.h>

#include <QLineF>
#include <QtMath>

#include "schrodinger/sketcher/rdkit/coord_utils.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"

#include "schrodinger/sketcher/molviewer/constants.h"

namespace schrodinger
{
namespace sketcher
{

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

QPointF best_placing_around_origin(const std::vector<QPointF>& points,
                                   const bool limit_to_120_for_single_neighbor)
{
    if (points.empty()) {
        return QPointF(BOND_LENGTH * VIEW_SCALE, 0);
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

    /**
     * find the biggest angle interval and return a point in the middle of it.
     * It is common that two or more intervals are nearly identical in size.
     * In that case, return the first one to avoid numerical instability in the
     * angle calculations.
     */
    const float TOLERANCE = 1.01;
    int best_i = 0;
    auto best_angle = (angles[best_i + 1] - angles[best_i]) * 0.5;
    for (unsigned int i = 0; i < angles.size() - 1; ++i) {
        auto angle_i = (angles[i + 1] - angles[i]) * 0.5;
        if (angle_i > best_angle * TOLERANCE) {
            best_i = i;
            best_angle = angle_i;
        }
    }
    // For a single substituent, limit the angle to 120 (instead of 180)
    bool limit_to_120 =
        limit_to_120_for_single_neighbor && (points.size() == 1);
    float max_angle = limit_to_120 ? 120.0 : 180.0;
    float angle_to_use = (angles[best_i + 1] - angles[best_i]) * 0.5;
    float angle_increment = std::min(max_angle, angle_to_use);
    float return_angle = angles[best_i] + angle_increment;
    QLineF line(origin, QPointF(BOND_LENGTH * VIEW_SCALE, 0));
    line.setAngle(return_angle);
    return line.p2();
}

RDGeom::Point3D
best_placing_around_origin(const std::vector<RDGeom::Point3D>& points,
                           const bool limit_to_120_for_single_neighbor)
{
    std::vector<QPointF> qpoints;
    qpoints.reserve(points.size());
    std::transform(points.begin(), points.end(), std::back_inserter(qpoints),
                   to_scene_xy);
    QPointF best =
        best_placing_around_origin(qpoints, limit_to_120_for_single_neighbor);
    return to_mol_xy(best);
}

RDGeom::Point3D flip_point(const RDGeom::Point3D& point,
                           const RDGeom::Point3D& start,
                           const RDGeom::Point3D& end)
{
    auto angle = get_angle_radians(end, start, point);
    if (std::isnan(angle)) {
        return point;
    }
    return rotate_point_radians(point, start, -2.0 * angle);
}

RDGeom::Point3D rotate_point(const RDGeom::Point3D& point,
                             const RDGeom::Point3D& center_of_rotation,
                             float angle)
{
    auto angle_radians = qDegreesToRadians(angle);
    return rotate_point_radians(point, center_of_rotation, angle_radians);
}

RDGeom::Point3D rotate_point_radians(const RDGeom::Point3D& point,
                                     const RDGeom::Point3D& center_of_rotation,
                                     float angle_radians)
{
    auto cosine = std::cos(angle_radians);
    auto sine = std::sin(angle_radians);

    auto coord = point - center_of_rotation;
    auto rotated_coord = coord;
    rotated_coord.x = coord.x * cosine - coord.y * sine;
    rotated_coord.y = coord.x * sine + coord.y * cosine;
    return rotated_coord + center_of_rotation;
}

RDGeom::Point3D find_centroid(
    const RDKit::ROMol& mol,
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)

{
    return find_centroid(mol.getConformer(), atoms, non_molecular_objects);
}

RDGeom::Point3D find_centroid(
    const RDKit::Conformer& conf,
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{
    std::vector<RDGeom::Point3D> positions;
    if (atoms.empty() && non_molecular_objects.empty()) {
        positions = conf.getPositions();
    } else {
        positions.reserve(atoms.size() + non_molecular_objects.size());
        for (const auto& atom : atoms) {
            positions.push_back(conf.getAtomPos(atom->getIdx()));
        }
        for (const auto* non_mol_obj : non_molecular_objects) {
            positions.push_back(non_mol_obj->getCoords());
        }
    }
    return rdkit_extensions::compute_centroid(positions);
}

RDGeom::Point3D find_centroid(
    const RDKit::ROMol& mol,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{
    auto all_atoms = mol.atoms();
    std::unordered_set<const RDKit::Atom*> all_atoms_set(all_atoms.begin(),
                                                         all_atoms.end());
    return find_centroid(mol, all_atoms_set, non_molecular_objects);
}

void center_mol_on(RDKit::ROMol& mol, const RDGeom::Point3D& center)
{
    if (!mol.getNumConformers()) {
        return;
    }
    auto centroid = find_centroid(mol);
    auto& conf = mol.getConformer();
    for (auto atom : mol.atoms()) {
        auto current_pos = conf.getAtomPos(atom->getIdx());
        conf.setAtomPos(atom->getIdx(), current_pos - centroid + center);
    }
}
void center_on_origin(RDKit::ROMol& mol)
{
    center_mol_on(mol, RDGeom::Point3D(0.0, 0.0, 0.0));
}

QPainterPath get_wavy_line_path(const int number_of_waves,
                                const qreal width_per_wave, const qreal height,
                                const qreal angle)
{
    qreal half_width = width_per_wave / 2.0;
    qreal start_x = -number_of_waves * half_width;
    QRectF half_wiggle_rect(start_x, -height / 2.0, half_width, height);

    QPainterPath path;
    path.moveTo(start_x, 0);
    for (int i = 0; i < number_of_waves; ++i) {
        path.arcTo(half_wiggle_rect, -180, -180);
        half_wiggle_rect.translate(half_width, 0);
        path.arcTo(half_wiggle_rect, -180, 180);
        half_wiggle_rect.translate(half_width, 0);
    }

    QTransform transform;
    transform.rotate(angle);
    return transform.map(path);
}

qreal get_attachment_point_line_angle(const RDKit::Atom* const ap_atom)
{
    auto& mol = ap_atom->getOwningMol();
    // attachment point atoms have exactly one neighbor
    RDKit::Atom* neighbor = *(mol.atomNeighbors(ap_atom).begin());
    auto& conf = mol.getConformer();
    auto ap_pos = conf.getAtomPos(ap_atom->getIdx());
    auto neighbor_pos = conf.getAtomPos(neighbor->getIdx());
    QLineF bond_line(to_scene_xy(neighbor_pos), to_scene_xy(ap_pos));
    return -bond_line.normalVector().angle();
}

double get_angle_radians(const RDGeom::Point3D& a, const RDGeom::Point3D& b,
                         const RDGeom::Point3D& c)
{
    auto vec_ba = a - b;
    auto vec_bc = c - b;
    return vec_ba.signedAngleTo(vec_bc);
}

void rotate_conformer_radians(const double angle_radians,
                              const RDGeom::Point3D& center_of_rotation,
                              RDKit::Conformer& conf)
{
    for (auto& cur_atom_coords : conf.getPositions()) {
        cur_atom_coords = rotate_point_radians(
            cur_atom_coords, center_of_rotation, angle_radians);
    }
}

bool are_points_on_same_side_of_line(const RDGeom::Point3D& point1,
                                     const RDGeom::Point3D& point2,
                                     const RDGeom::Point3D& line_endpoint)
{
    return are_points_on_same_side_of_line(
        to_scene_xy(point1), to_scene_xy(point2), to_scene_xy(line_endpoint));
}

bool are_points_on_same_side_of_line(const QPointF& point1,
                                     const QPointF& point2,
                                     const QPointF& line_endpoint)
{
    qreal x = line_endpoint.x();
    qreal y = line_endpoint.y();

    qreal slope, d1, d2;
    if (qFabs(x) > qFabs(y)) {
        slope = y / x;
        d1 = point1.y() - slope * point1.x();
        d2 = point2.y() - slope * point2.x();
    } else {
        // the line might be vertical (or close to vertical), so we flip the
        // axes to avoid divide by zero (or arithmetic underflow) in our
        // calculations
        slope = x / y;
        d1 = point1.x() - slope * point1.y();
        d2 = point2.x() - slope * point2.y();
    }

    // the sign of d1 and d2 tell which side of the line they're on
    return std::signbit(d1) == std::signbit(d2);
}

qreal get_rounded_angle_radians(const QPointF& start, const QPointF& end)
{
    qreal angle = QLineF(start, end).angle();
    angle = qDegreesToRadians(angle);
    int rounded = std::round(angle * 6.0 / M_PI);
    return rounded / 6.0 * M_PI;
}

void trim_line_to_rect(QLineF& line, const QRectF& rect, qreal min_length)
{
    // expand rect by 4 pixels in all directions
    const QMarginsF margins(ATOM_LABEL_MARGIN, ATOM_LABEL_MARGIN,
                            ATOM_LABEL_MARGIN, ATOM_LABEL_MARGIN);
    const QRectF enlarged_rect = rect + margins;

    QPointF inter_point; // the intersection point

    // if the line is completely inside the rect, skip the trimming
    if (enlarged_rect.contains(line.p1()) &&
        enlarged_rect.contains(line.p2())) {
        return;
    }
    if (enlarged_rect.contains(line.p1()) &&
        intersection_of_line_and_rect(line, enlarged_rect, inter_point)) {
        if (QLineF(inter_point, line.p2()).length() < min_length) {
            // trim to exactly min_length. We use a temporary line in the
            // opposite direction because setLength() moves p2 and we want to
            // move p1 here
            auto trimmed_line = QLineF(line.p2(), inter_point);
            trimmed_line.setLength(min_length);
            line.setP1(trimmed_line.p2());
        } else {
            line.setP1(inter_point);
        }
    } else if (enlarged_rect.contains(line.p2()) &&
               intersection_of_line_and_rect(line, enlarged_rect,
                                             inter_point)) {
        if (QLineF(line.p1(), inter_point).length() < min_length) {
            auto trimmed_line = QLineF(line.p1(), inter_point);
            trimmed_line.setLength(min_length);
            line.setP2(trimmed_line.p2());
        } else {
            line.setP2(inter_point);
        }
    }
}

bool intersection_of_line_and_rect(const QLineF& line, const QRectF& rect,
                                   QPointF& inter_point)
{
    // break the bounding rect up into four lines
    QLineF top(rect.topLeft(), rect.topRight());
    QLineF left(rect.topLeft(), rect.bottomLeft());
    QLineF bottom(rect.bottomLeft(), rect.bottomRight());
    QLineF right(rect.topRight(), rect.bottomRight());

    bool found_intersection =
        (line.intersects(top, &inter_point) == QLineF::BoundedIntersection ||
         line.intersects(left, &inter_point) == QLineF::BoundedIntersection ||
         line.intersects(bottom, &inter_point) == QLineF::BoundedIntersection ||
         line.intersects(right, &inter_point) == QLineF::BoundedIntersection);
    return found_intersection;
}

RDGeom::Point3D
move_molecule_to_the_right_of(RDKit::ROMol& mol_to_move,
                              const RDKit::ROMol& stationary_mol,
                              const double gap_size)
{
    auto& to_move_conf = mol_to_move.getConformer();
    const auto& stationary_conf = stationary_mol.getConformer();
    if (!stationary_conf.getNumAtoms() || !to_move_conf.getNumAtoms()) {
        throw std::runtime_error("No atoms found");
    }

    double max_stationary_x = -DBL_MAX;
    double center_stationary_y = 0.0;
    for (const auto& cur_coords : stationary_conf.getPositions()) {
        max_stationary_x = std::max(max_stationary_x, cur_coords.x);
        center_stationary_y += cur_coords.y;
    }
    center_stationary_y /= stationary_conf.getNumAtoms();

    double min_to_move_x = DBL_MAX;
    double center_to_move_y = 0.0;
    for (const auto& cur_coords : to_move_conf.getPositions()) {
        min_to_move_x = std::min(min_to_move_x, cur_coords.x);
        center_to_move_y += cur_coords.y;
    }
    center_to_move_y /= to_move_conf.getNumAtoms();

    const double delta_x = max_stationary_x - min_to_move_x + gap_size;
    const double delta_y = center_stationary_y - center_to_move_y;
    const RDGeom::Point3D offset(delta_x, delta_y, 0);
    for (auto& cur_coords : to_move_conf.getPositions()) {
        cur_coords += offset;
    }

    return RDGeom::Point3D(max_stationary_x + gap_size / 2, center_stationary_y,
                           0);
}

} // namespace sketcher
} // namespace schrodinger
