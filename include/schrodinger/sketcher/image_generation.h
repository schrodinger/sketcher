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

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/image_constants.h"

template <typename T> class QList;

namespace RDKit
{
class ROMol;
class ChemicalReaction;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * Supported image formats
 */
enum class ImageFormat { PNG, SVG };

/**
 * @param format file format to consider
 * @return supported image extension of the given format
 * @throw std::out_of_range if format is not supported
 */
SKETCHER_API std::string get_image_extension(ImageFormat format);

/**
 * Image generation options
 */

const qreal AUTOSCALE = -1.0;

struct RenderOptions {
    QSize width_height = QSize(400, 400);
    QColor background_color = Qt::transparent;
    // The scale used to size the molecule to width_height.  If < 0, the
    // molecule will be scaled automatically to fill the space.  If the scale is
    // too large for the current molecule (i.e if it would make the image larger
    // than width_height), then it will be ignored and auto scaling will be used
    // instead.
    qreal scale = AUTOSCALE;

    // If false, the generated image will be exactly width_height, but may
    // contain padding. (E.g., if width_height is 200x1000 for a square-ish
    // molecule, then the generated image will contain large amounts of
    // horizontal padding.) If trim_image true, the generated image will be no
    // larger than width_height and will not contain any padding. (E.g., if
    // width_height is 200x1000 for a square-ish molecule, then the generated
    // image will be 200x~200.)
    bool trim_image = false;

    // The font size to use for atom labels.  Other fonts will be scaled
    // relative to this.
    int font_size = DEFAULT_FONT_SIZE;

    // Scaling to apply to bond width.  Values less than one will result in
    // thinner bonds, while values greater than one will result in thicker
    // bonds.
    qreal bond_width_scale = 1.0;

    // User annotations and colorings based on atom and bond indices.
    // Bonds between haloed atoms will not be haloed unless separately
    // specified. Bonds lines between atoms will have split color based
    // on the line color of connected atoms, unless bond line coloring is
    // separately specified.
    QHash<int, std::string> rdatom_index_to_label;
    QHash<int, QColor> rdatom_index_to_halo_color;
    QHash<int, QColor> rdbond_index_to_halo_color;
    QHash<int, QColor> rdatom_index_to_line_color;
    QHash<int, QColor> rdbond_index_to_line_color;
    // Molecular rendering options
    StereoLabels show_stereo_annotations = StereoLabels::ALL;
    bool show_absolute_stereo_groups = false;
    bool show_simplified_stereo_annotation = false;
    bool show_symbol_for_H_isotopes = false;

    // Note that the following settings are only fully supported when using
    // molviewer
    CarbonLabels carbon_labels = CarbonLabels::NONE;
    ColorScheme color_scheme = ColorScheme::DEFAULT;
};

/**
 * @param mol/rxn/text molecule to render
 * @param opts given image generation configuration
 * @return image generated from the 2D sketcher
 */
SKETCHER_API QPicture get_qpicture(const RDKit::ROMol& mol,
                                   const RenderOptions& opts = RenderOptions());
SKETCHER_API QPicture get_qpicture(const RDKit::ChemicalReaction& rxn,
                                   const RenderOptions& opts = RenderOptions());
SKETCHER_API QPicture get_qpicture(const std::string& text,
                                   const RenderOptions& opts = RenderOptions());

/**
 * @param mol/rxn/text molecule to render
 * @param opts given image generation configuration
 * @return picture generated from the 2D sketcher
 */
SKETCHER_API QImage get_qimage(const RDKit::ROMol& mol,
                               const RenderOptions& opts = RenderOptions());
SKETCHER_API QImage get_qimage(const RDKit::ChemicalReaction& rxn,
                               const RenderOptions& opts = RenderOptions());
SKETCHER_API QImage get_qimage(const std::string& text,
                               const RenderOptions& opts = RenderOptions());

/**
 * @param mol/rxn/text molecule to render
 * @param format format of the image
 * @param opts given image generation configuration
 * @return byte array of data generated from the 2D sketcher
 */
SKETCHER_API QByteArray
get_image_bytes(const RDKit::ROMol& mol, ImageFormat format,
                const RenderOptions& opts = RenderOptions());
SKETCHER_API QByteArray
get_image_bytes(const RDKit::ChemicalReaction& rxn, ImageFormat format,
                const RenderOptions& opts = RenderOptions());
SKETCHER_API QByteArray
get_image_bytes(const std::string& text, ImageFormat format,
                const RenderOptions& opts = RenderOptions());

// FIXME: Remove this once image_generation is fully switched over to molviewer
SKETCHER_API QByteArray
get_LiveDesign_image_bytes(const RDKit::ROMol& mol, ImageFormat format,
                           const RenderOptions& opts = RenderOptions());
SKETCHER_API QByteArray get_LiveDesign_image_bytes(
    const RDKit::ChemicalReaction& rxn, ImageFormat format,
    const RenderOptions& opts = RenderOptions());

/**
 * @param filepath path to stock image to generate should it be known that the
 * 2D sketcher image either can't or shouldn't render the image. Used primarily
 * for LiveDesign to create entity images for large format structure classes.
 */
SKETCHER_API QByteArray
get_stock_image_bytes(const std::string& filepath, ImageFormat format,
                      const RenderOptions& opts = RenderOptions());

/**
 * @param mol/rxn/text molecule to render
 * @param filename path to write to; format interpreted from extension
 * @param opts given image generation configuration
 */
SKETCHER_API void save_image_file(const RDKit::ROMol& mol,
                                  const std::string& filename,
                                  const RenderOptions& opts = RenderOptions());
SKETCHER_API void save_image_file(const RDKit::ChemicalReaction& rxn,
                                  const std::string& filename,
                                  const RenderOptions& opts = RenderOptions());
SKETCHER_API void save_image_file(const std::string& text,
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

/**
 * Get the scale that will be used for rendering the provided mol when
 * RenderOptions::scale is set to AUTOSCALE.
 * @param rdmol The molecule to be rendered
 * @param opts The image generation options to use.  Note that opts.scale will
 * be ignored.
 * @return The scaling value
 */
SKETCHER_API qreal get_autoscale_value(
    const RDKit::ROMol& rdmol, const RenderOptions& opts = RenderOptions());

/**
 * Get the scale that will be used for rendering the provided reaction when
 * RenderOptions::scale is set to AUTOSCALE.
 * @param rxn The reaction to be rendered
 * @param opts The image generation options to use.  Note that opts.scale will
 * be ignored.
 * @return The scaling value
 */
SKETCHER_API qreal
get_autoscale_value(const RDKit::ChemicalReaction& rxn,
                    const RenderOptions& opts = RenderOptions());

} // namespace sketcher
} // namespace schrodinger
