/* -------------------------------------------------------------------------
 * Implements schrodinger::sketcher:: rendering APIs
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include "schrodinger/sketcher/image_generation.h"

#include <QBuffer>
#include <QByteArray>
#include <QFile>
#include <QFileInfo>
#include <QPainter>
#include <boost/algorithm/string/predicate.hpp>
#include <qsvggenerator.h>

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
    for (auto index_color : opts.highlight_atom_index_to_color) {
        atoms.at(index_color.first)->setHighlightColor(index_color.second);
    }
    auto bonds = scene.quickGetBonds();
    for (auto bond : bonds) {
        bond->removeHighlight();
    }
    for (auto index_color : opts.highlight_bond_index_to_color) {
        bonds.at(index_color.first)->setHighlightColor(index_color.second);
    }
}

void paint_scene(QPaintDevice* device, const RDKit::ROMol& rdmol,
                 const RenderOptions& opts)
{
    sketcherScene scene;
    SketcherModel sketcher_model(&scene);
    scene.setModel(&sketcher_model);

    scene.addRDKitMolecule(rdmol);
    setHighlights(scene, opts);

    scene._backgroundColor = Qt::transparent;
    auto scene_rect = scene.itemsBoundingRect();

    QPainter painter(device);
    auto target_rect = painter.viewport();
    painter.eraseRect(target_rect);
    painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    painter.fillRect(target_rect, opts.background_color);
    scene.render(&painter, target_rect, scene_rect);
}

ImageFormat get_format(const QString& filename)
{
    const std::unordered_map<std::string, ImageFormat> ext_to_format_map = {
        {"png", ImageFormat::PNG}, {"svg", ImageFormat::SVG}};
    auto ext = QFileInfo(filename).suffix().toStdString();
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
    paint_scene(&picture, rdmol, opts);
    return picture;
}

QImage get_qimage(const RDKit::ROMol& rdmol, const RenderOptions& opts)
{
    QImage image(opts.width_height, QImage::Format_ARGB32);
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

} // namespace sketcher
} // namespace schrodinger
