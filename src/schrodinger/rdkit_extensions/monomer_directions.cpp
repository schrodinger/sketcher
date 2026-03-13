#include "schrodinger/rdkit_extensions/monomer_directions.h"

#include <map>
#include <numbers>
#include <stdexcept>

namespace schrodinger::rdkit_extensions
{

constexpr double SQRT1_2 = 1.0 / std::numbers::sqrt2;

static const std::map<Direction, RDGeom::Point3D> DIRECTION_TO_UNIT_VECTOR_MAP =
    {{Direction::N, RDGeom::Point3D(0, 1, 0)},
     {Direction::S, RDGeom::Point3D(0, -1, 0)},
     {Direction::E, RDGeom::Point3D(1, 0, 0)},
     {Direction::W, RDGeom::Point3D(-1, 0, 0)},
     {Direction::NE, RDGeom::Point3D(SQRT1_2, SQRT1_2, 0)},
     {Direction::NW, RDGeom::Point3D(-SQRT1_2, SQRT1_2, 0)},
     {Direction::SE, RDGeom::Point3D(SQRT1_2, -SQRT1_2, 0)},
     {Direction::SW, RDGeom::Point3D(-SQRT1_2, -SQRT1_2, 0)}};

RDGeom::Point3D direction_to_unit_vector(Direction dir)
{
    auto it = DIRECTION_TO_UNIT_VECTOR_MAP.find(dir);
    if (it != DIRECTION_TO_UNIT_VECTOR_MAP.end()) {
        return it->second;
    }
    throw std::runtime_error("Invalid direction");
}

} // namespace schrodinger::rdkit_extensions
