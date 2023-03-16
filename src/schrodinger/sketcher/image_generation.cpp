/* -------------------------------------------------------------------------
 * Implements schrodinger::sketcher:: rendering APIs
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
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

#include "schrodinger/sketcher/qt_utils.h"
#include "schrodinger/sketcher/Scene.h"
#include "schrodinger/sketcher/sketcher.h"
#include "schrodinger/sketcher/sketcher_model.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

void setHighlights(const sketcherScene& scene, const RenderOptions& opts)
{
    auto atoms = scene.quickGetAtoms();
    for (auto atom : atoms) {
        atom->removeHighlight();
    }
    for (auto [index, color] :
         asKeyValue(opts.rdatom_index_to_highlight_color)) {
        atoms.at(index)->setHighlightColor(color);
    }
    auto bonds = scene.quickGetBonds();
    for (auto bond : bonds) {
        bond->removeHighlight();
    }
    for (auto [index, color] :
         asKeyValue(opts.rdbond_index_to_highlight_color)) {
        bonds.at(index)->setHighlightColor(color);
    }
}

void setUserAnnotations(const sketcherScene& scene, const RenderOptions& opts)
{
    auto atoms = scene.quickGetAtoms();
    for (auto atom : atoms) {
        atom->removeUserAnnotation();
    }
    for (auto [index, text] : asKeyValue(opts.rdatom_index_to_annotation)) {
        atoms.at(index)->setUserAnnotation(text);
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

void paint_scene(QPaintDevice* device, const RDKit::ROMol& rdmol,
                 const RenderOptions& opts)
{
    sketcherScene scene;
    SketcherModel sketcher_model(&scene);
    scene.setModel(&sketcher_model);

    scene.addRDKitMolecule(rdmol);
    setHighlights(scene, opts);
    setUserAnnotations(scene, opts);

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

} // unnamed namespace

QPicture get_qpicture(const RDKit::ROMol& rdmol, const RenderOptions& opts)
{
    QPicture picture;
    picture.setBoundingRect(QRect(QPoint(0, 0), opts.width_height));
    paint_scene(&picture, rdmol, opts);
    return picture;
}

QImage get_qimage(const RDKit::ROMol& rdmol, const RenderOptions& opts)
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
    paint_scene(&image, rdmol, opts);
    return image;
}

QByteArray get_image_bytes(const RDKit::ROMol& rdmol, ImageFormat format,
                           const RenderOptions& opts)
{
    QBuffer buffer;
    buffer.open(QIODevice::WriteOnly);

    if (format == ImageFormat::PNG) {
        auto image = get_qimage(rdmol, opts);
        image.save(&buffer, "PNG");

    } else if (format == ImageFormat::SVG) {
        QSvgGenerator svg_gen;
        svg_gen.setSize(opts.width_height);
        svg_gen.setViewBox(QRect(QPoint(0, 0), opts.width_height));
        svg_gen.setOutputDevice(&buffer);
        paint_scene(&svg_gen, rdmol, opts);

    } else {
        throw std::runtime_error("Unknown ImageFormat");
    }

    return buffer.data();
}

void save_image_file(const RDKit::ROMol& rdmol, const std::string& filename,
                     const RenderOptions& opts)
{
    auto path = QString::fromStdString(filename);
    auto format = get_format(path);
    auto data = get_image_bytes(rdmol, format, opts);
    QFile file(path);
    file.open(QIODevice::WriteOnly);
    file.write(data);
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
