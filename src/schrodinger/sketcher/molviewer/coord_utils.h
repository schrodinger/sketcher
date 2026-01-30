#pragma once

#include <vector>
#include <unordered_set>

#include <QtGlobal>
#include <QPainterPath>
#include <QPointF>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"

namespace RDGeom
{
class Point3D;
} // namespace RDGeom

namespace RDKit
{
class Atom;
class Conformer;
class ROMol;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * If the given molecule has coordinates, center it on the origin. Otherwise do
 * nothing.
 */
void center_mol_on(RDKit::ROMol& mol, const RDGeom::Point3D& center);

void center_on_origin(RDKit::ROMol& mol);

/**
 * Covert a location from MolModel/RDKit coordinates to Scene coordinates.
 * @param xyz RDKit atom coordinates to transform
 * @return resulting 2D scene coordinates
 */
SKETCHER_API QPointF to_scene_xy(const RDGeom::Point3D& mol_xy);

/**
 * Covert a location from Scene coordinates to MolModel/RDKit coordinates.
 * @param scene_xy Scene coordinates to transform
 * @return resulting MolModel coordinates
 */
SKETCHER_API RDGeom::Point3D to_mol_xy(const QPointF& scene_xy);

/**
 * For each neighbor of the given atom, return the coordinates of the neighbor
 * relative to the given atom (i.e. assume that atom is at (0, 0)).
 * @return The resulting relative MolModel coordinates
 */
std::vector<RDGeom::Point3D>
get_relative_positions_of_atom_neighbors(const RDKit::Atom* const atom);

/**
 * Determine the best location to place something new around an atom.
 * @param points All points to take into account, given in Scene coordinates
 * relative to the atom.  These are typically coordinates of neighboring atoms
 * or labels.
 * @param limit_to_120_for_single_neighbor If points contains only a single
 * point, should we limit the angle between this point and the returned point to
 * be 120 degrees or less.  If false, the angle can be up to 180 degrees.
 * @return A point one bond length (i.e. one unit in the RDKit coordinate
 * system) away from the origin, in the direction that places it furthest from
 * any of the given points.  Given in Scene coordinates relative to the atom.
 */
SKETCHER_API QPointF
best_placing_around_origin(const std::vector<QPointF>& points,
                           const bool limit_to_120_for_single_neighbor = true);

/**
 * An overload of the above method that operates in relative MolModel
 * coordinates instead of relative Scene coordinates.
 */
SKETCHER_API RDGeom::Point3D
best_placing_around_origin(const std::vector<RDGeom::Point3D>& points,
                           const bool limit_to_120_for_single_neighbor = true);

/**
 * Create and return a path for painting a wavy line
 *
 * @param number_of_waves The number of "waves" in the squiggle, where each wave
 * goes from zero, to the maximum, to the minimum, and back to
 * zero
 * @param width_per_wave The width of each "wave," measured in Scene units.
 * @param height The height of the wavy line, measured in Scene units
 * @param angle The angle of wavy line, measured in degrees
 */
SKETCHER_API QPainterPath get_wavy_line_path(const int number_of_waves,
                                             const qreal width_per_wave,
                                             const qreal height,
                                             const qreal angle);

/**
 * Determine the correct angle for painting a squiggle used to represent an
 * attachment point
 *
 * @param atom The attachment point atom
 * @return The angle, measured in degrees
 */
SKETCHER_API qreal
get_attachment_point_line_angle(const RDKit::Atom* const atom);

/**
 * @param point coordinate to rotate
 * @param center_of_rotation center about which to rotate
 * @param angle specified rotation angle in degrees
 * @return a new point generated from the requested rotation
 */
SKETCHER_API RDGeom::Point3D
rotate_point(const RDGeom::Point3D& point,
             const RDGeom::Point3D& center_of_rotation, float angle);

/**
 * A version of rotate_point that takes the angle in radians instead of degrees
 */
SKETCHER_API RDGeom::Point3D
rotate_point_radians(const RDGeom::Point3D& point,
                     const RDGeom::Point3D& center_of_rotation,
                     float angle_radians);

/**
 * mirror coordinates about a segment
 * @param point coordinate to flip
 * @param start start of the segment about which to flip
 * @param end end of the segment about which to flip
 */
SKETCHER_API RDGeom::Point3D flip_point(const RDGeom::Point3D& point,
                                        const RDGeom::Point3D& start,
                                        const RDGeom::Point3D& end);

/**
 * @return the centroid of a set of atoms and non-molecular objects. If no
 * atoms or non-molecular objects are given, the centroid of the whole
 * molecule plus all non-molecular objects is returned.
 * @param mol the molecule to compute the centroid for
 * @param atoms the atoms to compute the centroid for
 * @param non_molecular_objects the non-molecular objects to consider in
 * addition to mol/atoms
 */
SKETCHER_API RDGeom::Point3D
find_centroid(const RDKit::ROMol& mol,
              const std::unordered_set<const RDKit::Atom*>& atoms = {},
              const std::unordered_set<const NonMolecularObject*>&
                  non_molecular_objects = {});

/**
 * An overload of find_centroid that takes a conformer in place of a molecule
 *
 * @overload
 */
SKETCHER_API RDGeom::Point3D
find_centroid(const RDKit::Conformer& conf,
              const std::unordered_set<const RDKit::Atom*>& atoms = {},
              const std::unordered_set<const NonMolecularObject*>&
                  non_molecular_objects = {});

/**
 * @return the centroid of a set of all atoms in a molecule plus additional
 * non-molecular objects.
 * @param mol the molecule to compute the centroid for
 * @param non_molecular_objects the non-molecular objects
 *
 * @overload
 */
SKETCHER_API RDGeom::Point3D find_centroid(
    const RDKit::ROMol& mol,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects);

/**
 * @return the angle in radians between vectors ab and bc
 */
SKETCHER_API double get_angle_radians(const RDGeom::Point3D& a,
                                      const RDGeom::Point3D& b,
                                      const RDGeom::Point3D& c);

/**
 * Rotate a conformer by the specified angle
 * @param[in] angle_radians The angle in radians
 * @param[in] center_of_rotation The point to rotate about
 * @param[in,out] conf The conformer to rotate.  This conformer will be modified
 * in place
 */
SKETCHER_API void
rotate_conformer_radians(const double angle_radians,
                         const RDGeom::Point3D& center_of_rotation,
                         RDKit::Conformer& conf);

/**
 * Determine whether two points are on the same or different sides of
 * a line.
 * @param point1 The first point
 * @param point2 The second point
 * @param line_endpoint A point on the line.  The other end of the line is
 * assumed to be (0, 0).
 * @return true if the points are on the same side of the line.  false
 * otherwise.
 */
SKETCHER_API bool
are_points_on_same_side_of_line(const RDGeom::Point3D& point1,
                                const RDGeom::Point3D& point2,
                                const RDGeom::Point3D& line_endpoint);

/**
 * An overload that accepts QPointFs instead of RDGeom::Point3Ds
 * @overload
 */
SKETCHER_API bool are_points_on_same_side_of_line(const QPointF& point1,
                                                  const QPointF& point2,
                                                  const QPointF& line_endpoint);

/**
 * @return the angle in radians of a line from start to end, rounded to the
 * nearest pi/6 (30 degree) increment.
 */
SKETCHER_API qreal get_rounded_angle_radians(const QPointF& start,
                                             const QPointF& end);

/**
 * Trim a line so that it ends at least 4 pixels outside of the given
 * rectangle.  Note that this method assumes that the line is either
 * completely outside of the rectangle (in which case no trimming is
 * required), or that the line has exactly one endpoint within the rectangle
 * (in which case that endpoint is moved outside of the rectangle).  We
 * assume that the line is *not* completely contained within the rectangle
 * and that the line does not pass through the rectangle and come out the
 * other side.
 *
 * @param line[in,out] The line to trim
 * @param subrect[in] The rectangle to trim to
 * @param min_length[in] The minimum length the line should be.  The line will
 * not be trimmed to less than this.
 */
SKETCHER_API void trim_line_to_rect(QLineF& line, const QRectF& subrect,
                                    qreal min_length = 0);

/**
 * Calculate the intersection point (if any) of a line and a rectangle.  See
 * caveats in the `trim_line_to_rect` docstring.
 *
 * @param line[in] The line
 * @param rect[in] The rectangle
 * @param i[out] Will be set to the intersection point, if there is one
 * @return Whether the line intersects the rectangle
 */
SKETCHER_API bool intersection_of_line_and_rect(const QLineF& line,
                                                const QRectF& rect, QPointF& i);

/**
 * Update the coordinates of a molecule so that it's directly to the right of
 * another molecule.
 * @param[in,out] mol_to_move The molecule to move
 * @param[in] stationary_mol The molecule that mol_to_move is moved to the
 * right of
 * @param[in] gap_size The desired distance between the right side of
 * stationary_mol and the left side of mol_to_move
 * @return The coordinates of the midpoint between the right side of
 * stationary_mol and the left side of mol_to_move
 * @throw std::runtime_error if either mol_to_move or stationary_mol have no
 * atoms
 */
SKETCHER_API RDGeom::Point3D
move_molecule_to_the_right_of(RDKit::ROMol& mol_to_move,
                              const RDKit::ROMol& stationary_mol,
                              const double gap_size);

} // namespace sketcher
} // namespace schrodinger
