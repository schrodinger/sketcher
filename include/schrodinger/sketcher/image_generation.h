/* -------------------------------------------------------------------------
 * Declares schrodinger::sketcher:: rendering APIs
 --------------------------------------------------------------------------- */

#pragma once

#include <QByteArray>
#include <QColor>
#include <QHash>
#include <QImage>
#include <QPicture>
#include <QSize>
#include <unordered_map>

#include "schrodinger/sketcher/definitions.h"

template <typename T> class QList;

namespace RDKit
{
class ROMol;
}

namespace schrodinger
{
namespace sketcher
{

/**
 * Supported image formats
 */
enum class ImageFormat { PNG, SVG };

/**
 * Image generation options
 */

const qreal AUTOSCALE = -1.0;

struct RenderOptions {
    QSize width_height = QSize(400, 400);
    QColor background_color = Qt::transparent;
    /**
     * The scale used to size the molecule to width_height.  If < 0, the
     * molecule will be scaled automatically to fill the space.  If the scale is
     * too large for the current molecule (i.e if it would make the image larger
     * than width_height), then it will be ignored and auto scaling will be used
     * instead.
     */
    qreal scale = AUTOSCALE;
    QHash<int, std::string> rdatom_index_to_atom_label;
    QHash<int, QColor> rdatom_index_to_highlight_color;
    QHash<int, QColor> rdbond_index_to_highlight_color;
    // TODO: incorporate both RendererSettings and additional LiveDesign options
};

// TODO: Overload functions to allow std::string as input

/**
 * @param rdmol molecule to render
 * @param opts given image generation configuration
 * @return image generated from the 2D sketcher
 */
SKETCHER_API QPicture get_qpicture(const RDKit::ROMol& rdmol,
                                   const RenderOptions& opts = RenderOptions());

/**
 * @param rdmol molecule to render
 * @param opts given image generation configuration
 * @return picture generated from the 2D sketcher
 */
SKETCHER_API QImage get_qimage(const RDKit::ROMol& rdmol,
                               const RenderOptions& opts = RenderOptions());

/**
 * @param rdmol molecule to render
 * @param format format of the image
 * @param opts given image generation configuration
 * @return byte array of data generated from the 2D sketcher
 */
SKETCHER_API QByteArray
get_image_bytes(const RDKit::ROMol& rdmol, ImageFormat format,
                const RenderOptions& opts = RenderOptions());

/**
 * @param rdmol molecule to render
 * @param filename path to write to; format interpreted from extension
 * @param opts given image generation configuration
 */
SKETCHER_API void save_image_file(const RDKit::ROMol& rdmol,
                                  const std::string& filename,
                                  const RenderOptions& opts = RenderOptions());

/**
 * Get the correct scale for rendering all of the provided mols to identically
 * sized images.  To render the the molecules using this scaling, set
 * RenderOptions::scale to the return value.
 * @param all_rdmols All molecules to be rendered
 * @param opts The image generation options to use.  Note that opts.scale will
 * be ignored.
 * @return The maximum scale value that will allow all molecules in all_rdmols
 * to fit within an image of size opts.width_height
 */
SKETCHER_API qreal
get_best_image_scale(const QList<RDKit::ROMol*> all_rdmols,
                     const RenderOptions& opts = RenderOptions());

} // namespace sketcher
} // namespace schrodinger
