#pragma once

#include <vector>

#include <QtGlobal>
#include <QPainterPath>
#include <QPointF>

#include "schrodinger/sketcher/definitions.h"

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
 * @return whether a 2D Conformer is available for the given molecule.
 */
SKETCHER_API bool has_2d_conformer(const RDKit::ROMol& mol);

/**
 * @return the first 2D Conformer found on the given molecule.
 * If none is found, one is generated and returned.
 */
SKETCHER_API RDKit::Conformer& get_2d_conformer(RDKit::ROMol& mol);

/**
 * Rescale a conformer so that the most frequent bond length matches the given
 * reference bond length.
 */
SKETCHER_API void rescale(RDKit::Conformer& conformer, RDKit::ROMol& mol,
                          double reference_bond_length);

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
 * @return A point one bond length (i.e. one unit in the RDKit coordinate
 * system) away from the origin, in the direction that places it furthest from
 * any of the given points.  Given in Scene coordinates relative to the atom.
 */
SKETCHER_API QPointF
best_placing_around_origin(const std::vector<QPointF>& points);

/**
 * An overload of the above method that operates in relative MolModel
 * coordinates instead of relative Scene coordinates.
 */
SKETCHER_API RDGeom::Point3D
best_placing_around_origin(const std::vector<RDGeom::Point3D>& points);

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

} // namespace sketcher
} // namespace schrodinger
