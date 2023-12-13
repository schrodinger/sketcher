/* -------------------------------------------------------------------------
 * Implements schrodinger::sketcher:: rendering APIs
 --------------------------------------------------------------------------- */

#include "schrodinger/sketcher/image_generation.h"

#include <limits>

#include <QBuffer>
#include <QByteArray>
#include <QFile>
#include <QFileInfo>
#include <QList>
#include <QPainter>
#include <boost/algorithm/string/predicate.hpp>
#include <qsvggenerator.h>

#include "schrodinger/sketcher/Scene.h"
#include "schrodinger/sketcher/highlighting_item.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/sketcher.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * Wrapper for QHash/QMap to use range-based for loops and structured bindings
 */
template <typename T> class asKeyValue
{
  public:
    asKeyValue(const T& data) : m_data{data}
    {
    }
    const auto begin()
    {
        return m_data.keyValueBegin();
    }
    const auto end()
    {
        return m_data.keyValueEnd();
    }

  private:
    const T& m_data;
};

void setHighlights(const sketcherScene& scene, const RenderOptions& opts)
{
    auto atoms = scene.quickGetAtoms();
    for (auto atom : atoms) {
        atom->removeHighlight();
    }
    for (auto [index, color] : asKeyValue(opts.rdatom_index_to_halo_color)) {
        atoms.at(index)->setHighlightColor(color);
    }
    auto bonds = scene.quickGetBonds();
    for (auto bond : bonds) {
        bond->removeHighlight();
    }
    for (auto [index, color] : asKeyValue(opts.rdbond_index_to_halo_color)) {
        bonds.at(index)->setHighlightColor(color);
    }
}

void setUserLabels(const sketcherScene& scene, const RenderOptions& opts)
{
    auto atoms = scene.quickGetAtoms();
    for (auto atom : atoms) {
        atom->removeUserLabel();
    }
    for (auto [index, text] : asKeyValue(opts.rdatom_index_to_label)) {
        atoms.at(index)->setUserLabel(text);
    }
}

qreal get_scale(const QRectF& scene_rect, const QSize& render_size)
{
    qreal scene_width = scene_rect.width();
    qreal scene_height = scene_rect.height();
    qreal x_ratio = render_size.width() / scene_width;
    qreal y_ratio = render_size.height() / scene_height;
    return qMin(x_ratio, y_ratio);
}

/**
 * @internal
 * Helper functions for the templated image generation APIs to inject various
 * inputs into the given sketcher scene.
 */
void add_to_scene(sketcherScene* scene, const RDKit::ROMol& rdmol)
{
    scene->addRDKitMolecule(rdmol);
}

void add_to_scene(sketcherScene* scene, const RDKit::ChemicalReaction& rxn)
{
    scene->addRDKitReaction(rxn);
}

void add_to_scene(sketcherScene* scene, const std::string& text)
{
    scene->importText(text);
}

template <typename T> void paint_scene(QPaintDevice* device, const T& input,
                                       const RenderOptions& opts)
{
    sketcherScene scene;
    SketcherModel sketcher_model(&scene);
    scene.setModel(&sketcher_model);

    add_to_scene(&scene, input);
    setHighlights(scene, opts);
    setUserLabels(scene, opts);

    // The atom and bond items use a semi-transparent version of the background
    // color to fade out bonds in order to make labels more readable.
    // Using Qt's black-transparent causes these to turn into black smudges.
    scene._backgroundColor = QColor(255, 255, 255, 0);

    QPainter painter(device);
    auto target_rect = painter.viewport();
    painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    painter.fillRect(target_rect, opts.background_color);

    // center the scene within the painter's viewport
    auto scene_rect = scene.findBoundingRect();
    qreal scale = get_scale(scene_rect, opts.width_height);
    if (opts.scale > 0 && opts.scale < scale) {
        // if the user has specified a scale, use that unless it would make the
        // molecule too big for the image
        scale = opts.scale;
    }
    QRectF centered_rect(0, 0, scene_rect.width() * scale,
                         scene_rect.height() * scale);
    centered_rect.moveCenter(target_rect.center());

    scene.render(&painter, centered_rect, scene_rect);
}

ImageFormat get_format(const QString& filename)
{
    const std::unordered_map<std::string, ImageFormat> ext_to_format_map = {
        {"png", ImageFormat::PNG}, {"svg", ImageFormat::SVG}};
    auto ext = QFileInfo(filename).completeSuffix().toStdString();
    try {
        return ext_to_format_map.at(ext);
    } catch (const std::out_of_range&) {
        throw std::invalid_argument("Unsupported file extension: " +
                                    filename.toStdString());
    }
}

template <typename T>
QPicture get_qpicture(const T& input, const RenderOptions& opts)
{
    QPicture picture;
    picture.setBoundingRect(QRect(QPoint(0, 0), opts.width_height));
    paint_scene(&picture, input, opts);
    return picture;
}

template <typename T>
QImage get_qimage(const T& input, const RenderOptions& opts)
{
    QImage image(opts.width_height, QImage::Format_ARGB32);
    {
        // initialize the contents of image. paint_scene will paint the correct
        // background color, so here we just need to replace uninitialized data
        // with transparent pixels.
        QPainter painter(&image);
        QRect target_rect(QPoint(0, 0), opts.width_height);
        painter.setCompositionMode(QPainter::CompositionMode_Source);
        painter.fillRect(target_rect, Qt::transparent);
    }
    paint_scene(&image, input, opts);
    return image;
}

template <typename T> QByteArray
get_image_bytes(const T& input, ImageFormat format, const RenderOptions& opts)
{
    QBuffer buffer;
    buffer.open(QIODevice::WriteOnly);

    if (format == ImageFormat::PNG) {
        auto image = get_qimage(input, opts);
        image.save(&buffer, "PNG");

    } else if (format == ImageFormat::SVG) {
        QSvgGenerator svg_gen;
        svg_gen.setSize(opts.width_height);
        svg_gen.setViewBox(QRect(QPoint(0, 0), opts.width_height));
        svg_gen.setOutputDevice(&buffer);
        paint_scene(&svg_gen, input, opts);

    } else {
        throw std::runtime_error("Unknown ImageFormat");
    }

    return buffer.data();
}

template <typename T> void save_image_file(const T& input,
                                           const std::string& filename,
                                           const RenderOptions& opts)
{
    auto path = QString::fromStdString(filename);
    auto format = get_format(path);
    auto data = get_image_bytes(input, format, opts);
    QFile file(path);
    file.open(QIODevice::WriteOnly);
    file.write(data);
}

} // unnamed namespace

QPicture get_qpicture(const RDKit::ROMol& mol, const RenderOptions& opts)
{
    return get_qpicture<RDKit::ROMol>(mol, opts);
}

QPicture get_qpicture(const RDKit::ChemicalReaction& rxn,
                      const RenderOptions& opts)
{
    return get_qpicture<RDKit::ChemicalReaction>(rxn, opts);
}

QPicture get_qpicture(const std::string& text, const RenderOptions& opts)
{
    return get_qpicture<std::string>(text, opts);
}

QImage get_qimage(const RDKit::ROMol& mol, const RenderOptions& opts)
{
    return get_qimage<RDKit::ROMol>(mol, opts);
}

QImage get_qimage(const RDKit::ChemicalReaction& rxn, const RenderOptions& opts)
{
    return get_qimage<RDKit::ChemicalReaction>(rxn, opts);
}

QImage get_qimage(const std::string& text, const RenderOptions& opts)
{
    return get_qimage<std::string>(text, opts);
}

QByteArray get_image_bytes(const RDKit::ROMol& mol, ImageFormat format,
                           const RenderOptions& opts)
{
    return get_image_bytes<RDKit::ROMol>(mol, format, opts);
}

QByteArray get_image_bytes(const RDKit::ChemicalReaction& rxn,
                           ImageFormat format, const RenderOptions& opts)
{
    return get_image_bytes<RDKit::ChemicalReaction>(rxn, format, opts);
}

QByteArray get_image_bytes(const std::string& text, ImageFormat format,
                           const RenderOptions& opts)
{
    return get_image_bytes<std::string>(text, format, opts);
}

void save_image_file(const RDKit::ROMol& mol, const std::string& filename,
                     const RenderOptions& opts)
{
    save_image_file<RDKit::ROMol>(mol, filename, opts);
}

void save_image_file(const RDKit::ChemicalReaction& rxn,
                     const std::string& filename, const RenderOptions& opts)
{
    save_image_file<RDKit::ChemicalReaction>(rxn, filename, opts);
}

void save_image_file(const std::string& text, const std::string& filename,
                     const RenderOptions& opts)
{
    save_image_file<std::string>(text, filename, opts);
}

qreal get_best_image_scale(const QList<RDKit::ROMol*> all_rdmols,
                           const RenderOptions& opts)
{
    if (all_rdmols.empty()) {
        return AUTOSCALE;
    }
    qreal best_scale = std::numeric_limits<qreal>::max();
    for (auto rdmol : all_rdmols) {
        // We sometimes get graphical artifacts when clearing and reusing a
        // sketcherScene instance for multiple molecules.  To avoid that, we
        // create a new scene instance for each molecule.  We can presumably
        // switch this code to use scene.clearInteractiveItems() instead once we
        // switch to the the molviewer Scene class here.
        sketcherScene scene;
        SketcherModel sketcher_model(&scene);
        scene.setModel(&sketcher_model);
        scene.addRDKitMolecule(*rdmol);
        qreal cur_scale =
            get_scale(scene.findBoundingRect(), opts.width_height);
        if (cur_scale < best_scale) {
            best_scale = cur_scale;
        }
    }
    return best_scale;
}

} // namespace sketcher
} // namespace schrodinger
