#pragma once

#include <string>
#include <tuple>
#include <utility>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/fonts.h"
#include "schrodinger/sketcher/tool/standard_scene_tool_base.h"

class QPointF;
class QRectF;

namespace RDKit
{
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{

namespace sketcher
{

/**
 * @return the attachment point name to a QString after converting apostrophes
 * to Unicode primes.
 */
SKETCHER_API QString prep_attachment_point_name(const std::string& name);

/**
 * Position the given rectangle to label a monomer's attachment point
 * @param ap_label_rect The rectangle to position. It should already be sized
 * correctly for the attachment point label.
 * @param monomer_coords The coordinates of the monomer being labeled
 * @param bound_coords The coordinates of the other monomer involved in the bond
 */
SKETCHER_API void position_ap_label_rect(QRectF& ap_label_rect,
                                         const QPointF& monomer_coords,
                                         const QPointF& bound_coords);

/**
 * Create and return a label for monomer's attachment point where it's connected
 * to bound_monomer. The calling scope take ownership of the returned graphics
 * item.  Note that if ap_name is empty, no graphics item will be created and a
 * null pointer will be returned.
 */
SKETCHER_API QGraphicsItem* create_label_for_bound_attachment_point(
    const RDKit::Atom* const monomer, const RDKit::Atom* const bound_monomer,
    const bool is_secondary_connection, const std::string& ap_name,
    const QColor& color, const Fonts& fonts, const Scene* const scene);

/**
 * @overload accepts the graphics item for monomer in place of the scene
 */
SKETCHER_API QGraphicsItem* create_label_for_bound_attachment_point(
    const RDKit::Atom* const monomer, const RDKit::Atom* const bound_monomer,
    const bool is_secondary_connection, const std::string& ap_name,
    const QColor& color, const Fonts& fonts,
    const QGraphicsItem* const monomer_item);

/**
 * Create and return labels for the attachment points on both ends of the
 * specified connector. The calling scope takes ownership of all returned
 * graphics items. Note that only a single (centered) graphics item will be
 * returned for nucleic acid base pair connections.
 */
SKETCHER_API std::vector<QGraphicsItem*>
create_attachment_point_labels_for_connector(const RDKit::Bond* const connector,
                                             const bool is_secondary_connection,
                                             const QColor& color,
                                             const Fonts& fonts,
                                             const Scene* const scene);

/**
 * @overload accepts the graphics items for the monomers involved in the
 * connection in place of the scene
 */
SKETCHER_API std::vector<QGraphicsItem*>
create_attachment_point_labels_for_connector(
    const RDKit::Bond* const connector, const bool is_secondary_connection,
    const QColor& color, const Fonts& fonts,
    const QGraphicsItem* const begin_monomer_item,
    const QGraphicsItem* const end_monomer_item);

} // namespace sketcher
} // namespace schrodinger
