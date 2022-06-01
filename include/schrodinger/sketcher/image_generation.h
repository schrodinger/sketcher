/* -------------------------------------------------------------------------
 * Declares schrodinger::sketcher:: rendering APIs
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <QByteArray>
#include <QColor>
#include <QImage>
#include <QPicture>
#include <QSize>

#include "schrodinger/sketcher/definitions.h"

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
struct RenderOptions {
    QSize width_height = QSize(400, 400);
    QColor background_color = Qt::transparent;
    // SKETCH-942: std::unordered_map<QColor, std::vector<int>> highlight_atoms;
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

} // namespace sketcher
} // namespace schrodinger
