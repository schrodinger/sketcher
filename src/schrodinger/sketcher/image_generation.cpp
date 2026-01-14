/* -------------------------------------------------------------------------
 * Implements schrodinger::sketcher:: rendering APIs
 --------------------------------------------------------------------------- */

#include "schrodinger/sketcher/image_generation.h"

#include <limits>
#include <memory>
#include <functional>

#include <QBuffer>
#include <QByteArray>
#include <QFile>
#include <QFileInfo>
#include <QGraphicsSvgItem>
#include <QList>
#include <QPainter>
#include <QPixmap>
#include <QSvgGenerator>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/assign.hpp>
#include <boost/bimap.hpp>
#include <unordered_map>

#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/constants.h"

/// a function that returns a paint device at the specified size
using PaintDeviceFunction =
    std::function<std::shared_ptr<QPaintDevice>(const QSize&)>;

template <> struct std::hash<QColor> {
    std::size_t operator()(const QColor& c) const noexcept
    {
        return std::hash<unsigned int>{}(c.rgba());
    }
};
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

using ImageFormatBimapType = boost::bimap<ImageFormat, std::string>;

// clang-format off
const ImageFormatBimapType IMAGE_FORMAT_EXT_BIMAP =
        boost::assign::list_of<ImageFormatBimapType::relation>
    (ImageFormat::PNG, ".png")
    (ImageFormat::SVG, ".svg");
// clang-format on

/**
 * @internal
 * Helper function to inject line colors data from the given RenderOptions into
 * the given MolModel
 */
void setLineColors(MolModel& model, const RenderOptions& opts)
{
    auto mol = model.getMol();
    // clear all line colors
    for (auto atom : mol->atoms()) {
        atom->clearProp(USER_COLOR);
    }
    for (auto bond : mol->bonds()) {
        bond->clearProp(USER_COLOR);
    }
    // set line colors
    for (auto [index, color] : asKeyValue(opts.rdatom_index_to_line_color)) {
        mol->getAtomWithIdx(index)->setProp(USER_COLOR, color);
    }
    for (auto [index, color] : asKeyValue(opts.rdbond_index_to_line_color)) {
        mol->getBondWithIdx(index)->setProp(USER_COLOR, color);
    }
}

/**
 * @internal
 * Helper function to inject halo highlighting data from the given RenderOptions
 * into the given MolModel
 */
void setHaloHighlightings(MolModel& model, const RenderOptions& opts)
{
    model.clearHaloHighlighting();
    std::unordered_map<QColor,
                       std::pair<std::unordered_set<const RDKit::Atom*>,
                                 std::unordered_set<const RDKit::Bond*>>>
        color_to_atoms_and_bonds;
    auto mol = model.getMol();
    for (auto [index, color] : asKeyValue(opts.rdatom_index_to_halo_color)) {
        color_to_atoms_and_bonds[color].first.insert(
            mol->getAtomWithIdx(index));
    }
    for (auto [index, color] : asKeyValue(opts.rdbond_index_to_halo_color)) {
        color_to_atoms_and_bonds[color].second.insert(
            mol->getBondWithIdx(index));
    }
    for (const auto& [color, atoms_and_bonds] : color_to_atoms_and_bonds) {
        model.addHaloHighlighting(atoms_and_bonds.first, atoms_and_bonds.second,
                                  color);
    }
}

/**
 * @internal
 * Helper function to inject atom labels data from the given RenderOptions into
 * the given MolModel
 */
void setAtomLabels(MolModel& model, const RenderOptions& opts)
{
    auto mol = model.getMol();
    for (auto [index, text] : asKeyValue(opts.rdatom_index_to_label)) {
        auto atom = mol->getAtomWithIdx(index);
        atom->setProp(RDKit::common_properties::_displayLabel, text);
    }
}

ImageFormat get_format(const QString& filename)
{
    auto ext = QFileInfo(filename).completeSuffix().toStdString();
    return IMAGE_FORMAT_EXT_BIMAP.right.at("." + ext);
}

std::pair<qreal, qreal> get_x_and_y_scales(const QRectF& scene_rect,
                                           const QSize& render_size)
{
    qreal scene_width = scene_rect.width();
    qreal scene_height = scene_rect.height();
    qreal x_ratio = render_size.width() / scene_width;
    qreal y_ratio = render_size.height() / scene_height;
    return {x_ratio, y_ratio};
}

/**
 * @return the appropriate scaling factor to make scene_rect as large as
 * possible, but no larger than render_size
 */
qreal get_scale(const QRectF& scene_rect, const QSize& render_size)
{
    auto [x_ratio, y_ratio] = get_x_and_y_scales(scene_rect, render_size);
    return qMin(x_ratio, y_ratio);
}

/**
 * @return what render_size should be trimmed to if we're rendering scene_rect
 * into render_size and RenderOptions::trim_image is true.
 */
QSize get_trimmed_size(const QRectF& scene_rect, const QSize& render_size)
{
    auto [x_ratio, y_ratio] = get_x_and_y_scales(scene_rect, render_size);
    qreal trimmed_width, trimmed_height;
    if (x_ratio < y_ratio) {
        trimmed_width = render_size.width();
        trimmed_height = x_ratio * scene_rect.height();
    } else {
        trimmed_width = y_ratio * scene_rect.width();
        trimmed_height = render_size.height();
    }
    return {static_cast<int>(trimmed_width), static_cast<int>(trimmed_height)};
}

/**
 * @return the appropriate image size for rendering the given scene
 */
QSize get_image_size(const RenderOptions& opts, const QRectF& scene_rect)
{
    if (opts.trim_image) {
        return get_trimmed_size(scene_rect, opts.width_height);
    } else {
        return opts.width_height;
    }
}

/**
 * @internal
 * Helper functions for the templated image generation APIs to inject various
 * inputs into the given sketcher scene.
 */
void add_to_mol_model(MolModel& mol_model, const RDKit::ROMol& rdmol)
{
    mol_model.addMol(rdmol);
}

void add_to_mol_model(MolModel& mol_model, const RDKit::ChemicalReaction& rxn)
{
    mol_model.addReaction(rxn);
}

void add_to_mol_model(MolModel& mol_model, const std::string& text)
{
    add_text_to_mol_model(mol_model, text);
}

void paint_scene_to_given_paint_device(QPaintDevice* device,
                                       const QGraphicsScene& scene,
                                       const QRectF& scene_rect,
                                       const RenderOptions& opts)
{

    QPainter painter(device);
    auto target_rect = painter.viewport();
    painter.setRenderHints(QPainter::Antialiasing | QPainter::TextAntialiasing);
    painter.fillRect(target_rect, opts.background_color);

    // center the scene within the painter's viewport
    qreal scale = get_scale(scene_rect, opts.width_height);
    if (opts.scale > 0 && opts.scale < scale) {
        // if the user has specified a scale, use that unless it would make the
        // molecule too big for the image
        scale = opts.scale;
    }

    QRectF centered_rect(0, 0, scene_rect.width() * scale,
                         scene_rect.height() * scale);
    centered_rect.moveCenter(target_rect.center());

    // QGraphicsScene::render is non-const, but we have to jump through const
    // hoops to consistently pass a QGraphicsScene through get_image_bytes and
    // the various const mol/rxn/text interfaces. The former public API
    // specifies that the QGraphicsScene parameter is non-const, and the later
    // all create a local Scene object within the static image generation,
    // so removing the const-ness from this scene seems like an acceptable
    // compromise to avoid duping the templated code for a non-const Scene.
    const_cast<QGraphicsScene*>(&scene)->render(&painter, centered_rect,
                                                scene_rect);
}

template <typename T>
void init_molviewer_image(MolModel& mol_model, SketcherModel& sketcher_model,
                          const T& input, const RenderOptions& opts)
{
    add_to_mol_model(mol_model, input);
    setLineColors(mol_model, opts);
    setHaloHighlightings(mol_model, opts);
    setAtomLabels(mol_model, opts);
    sketcher_model.loadRenderOptions(opts);
}

/**
 * Construct a paint device and render the given input to it
 *
 * @param input What to render into the scene.  Can be an RDKit molecule, an
 * RDKit reaction, a string containing a molecule (such as a SMILES string or a
 * MOL block), a pre-populated Scene (see overload below), or a stock image (see
 * other overload below)
 * @param opts The settings for rendering the molecule
 * @param instantiate_paint_device_to_size A function that takes a QSize and
 * returns a QPaintDevice instance of the specified size
 */
template <typename T> std::shared_ptr<QPaintDevice>
paint_scene(const T& input, const RenderOptions& opts,
            const PaintDeviceFunction& instantiate_paint_device_to_size)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);
    SketcherModel sketcher_model;
    Scene scene(&mol_model, &sketcher_model);
    init_molviewer_image(mol_model, sketcher_model, input, opts);

    return paint_scene(scene, opts, instantiate_paint_device_to_size);
}

// Save the scene directly as it is; used to save an image from SketcherWidget
template <> std::shared_ptr<QPaintDevice>
paint_scene(const Scene& input, const RenderOptions& opts,
            const PaintDeviceFunction& instantiate_paint_device_to_size)
{
    auto scene_rect = input.getSceneItemsBoundingRect();
    scene_rect.adjust(-SAVED_PICTURE_PADDING, -SAVED_PICTURE_PADDING,
                      SAVED_PICTURE_PADDING, SAVED_PICTURE_PADDING);
    auto image_size = get_image_size(opts, scene_rect);
    auto paint_device = instantiate_paint_device_to_size(image_size);
    paint_scene_to_given_paint_device(paint_device.get(), input, scene_rect,
                                      opts);
    return paint_device;
}

// Inject a prepared image into the scene; used to generate stock images
template <> std::shared_ptr<QPaintDevice>
paint_scene(const QFileInfo& file_info, const RenderOptions& opts,
            const PaintDeviceFunction& instantiate_paint_device_to_size)
{
    QGraphicsScene scene;
    auto file_path = file_info.filePath();
    switch (get_format(file_path)) {
        case ImageFormat::PNG:
            scene.addPixmap(QPixmap(file_path));
        case ImageFormat::SVG: {
            auto svg_item = new QGraphicsSvgItem(file_path);
            scene.addItem(svg_item);
        }
    }
    // we never trim the image here since we're displaying a placeholder, not
    // an actual molecule
    auto paint_device = instantiate_paint_device_to_size(opts.width_height);
    paint_scene_to_given_paint_device(paint_device.get(), scene,
                                      scene.itemsBoundingRect(), opts);
    return paint_device;
}

std::shared_ptr<QPicture> instantiate_qpicture_to_size(const QSize& size)
{
    auto picture = std::make_shared<QPicture>();
    picture->setBoundingRect(QRect(QPoint(0, 0), size));
    return picture;
}

std::shared_ptr<QImage> instantiate_qimage_to_size(const QSize& size)
{
    auto image = std::make_shared<QImage>(size, QImage::Format_ARGB32);
    // paint_scene will paint the correct background color, so here we just need
    // to replace uninitialized data with transparent pixels.
    QPainter painter(image.get());
    QRect target_rect(QPoint(0, 0), size);
    painter.setCompositionMode(QPainter::CompositionMode_Source);
    painter.fillRect(target_rect, Qt::transparent);
    return image;
}

std::shared_ptr<QSvgGenerator>
instantiate_qsvggenerator_to_size(const QSize& size)
{
    auto svg_gen = std::make_shared<QSvgGenerator>();
    // this buffer is assigned to a smart pointer in get_image_bytes
    auto buffer = new QBuffer();
    buffer->open(QIODevice::WriteOnly);
    svg_gen->setOutputDevice(buffer);
    svg_gen->setSize(size);
    svg_gen->setViewBox(QRect(QPoint(0, 0), size));
    return svg_gen;
}

template <typename T>
QPicture get_qpicture(const T& input, const RenderOptions& opts)
{
    auto paint_device =
        paint_scene<T>(input, opts, instantiate_qpicture_to_size);
    return *dynamic_cast<QPicture*>(paint_device.get());
}

template <typename T>
QImage get_qimage(const T& input, const RenderOptions& opts)
{
    auto paint_device = paint_scene<T>(input, opts, instantiate_qimage_to_size);
    return *dynamic_cast<QImage*>(paint_device.get());
}

template <typename T> QByteArray
get_image_bytes(const T& input, ImageFormat format, const RenderOptions& opts)
{
    QBuffer buffer;

    if (format == ImageFormat::PNG) {
        auto image = get_qimage(input, opts);
        buffer.open(QIODevice::WriteOnly);
        image.save(&buffer, "PNG");
        buffer.close();

    } else if (format == ImageFormat::SVG) {
        auto paint_device =
            paint_scene<T>(input, opts, instantiate_qsvggenerator_to_size);
        auto svg_gen = dynamic_cast<QSvgGenerator*>(paint_device.get());
        // instantiate_qsvggenerator_to_size constructed a new QBuffer on the
        // heap during the paint_scene call.  We assign that QBuffer to a smart
        // pointer here so that it will get destroyed at the end of this
        // function.
        std::shared_ptr<QBuffer> svg_buffer;
        svg_buffer.reset(dynamic_cast<QBuffer*>(svg_gen->outputDevice()));
        auto svg_data = svg_buffer->data();

        // Qt defines svg size in mm, but our code needs it specified in pixels.
        // This is a hack to make sure we get the right size definition.
        auto start_of_size = svg_data.indexOf("svg width");
        auto end_of_size = svg_data.indexOf(" viewBox", start_of_size);
        if (start_of_size == -1 || end_of_size == -1) {
            throw std::runtime_error("Failed to find svg size definition");
        }
        svg_data.remove(start_of_size, end_of_size - start_of_size);
        auto height_in_pxls = QString("svg width=\"%1px\" height=\"%2px\"\n")
                                  .arg(opts.width_height.width())
                                  .arg(opts.width_height.height());

        svg_data.insert(start_of_size, height_in_pxls.toUtf8());

        // Remove Qt's default title and description elements for cleaner SVG
        QString svg_string = QString::fromUtf8(svg_data);
        svg_string.remove("<title>Qt SVG Document</title>\n");
        svg_string.remove("<desc>Generated with Qt</desc>\n");
        svg_data = svg_string.toUtf8();

        buffer.setData(svg_data);
        svg_buffer->close();

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
    bool success = file.open(QIODevice::WriteOnly);
    if (!success) {
        throw std::runtime_error("Could not open file for writing");
    }
    file.write(data);
}

template <typename T>
qreal get_image_scale_for_mol_or_rxn(const T& mol_or_rxn, RenderOptions opts)

{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);
    SketcherModel sketcher_model;
    Scene scene(&mol_model, &sketcher_model);
    opts.scale = AUTOSCALE;
    init_molviewer_image(mol_model, sketcher_model, mol_or_rxn, opts);
    return get_scale(scene.getSceneItemsBoundingRect(), opts.width_height);
}

} // unnamed namespace

std::string get_image_extension(ImageFormat format)
{
    return IMAGE_FORMAT_EXT_BIMAP.left.at(format);
}

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

/**
 * @internal
 * Note that the function is not declared in image_generation.h, but rather in
 * scene.h to avoid introducing a dependency on Scene in the public headers
 */
QByteArray get_image_bytes(Scene& scene, ImageFormat format,
                           const RenderOptions& opts)
{
    return get_image_bytes<Scene>(scene, format, opts);
}

QByteArray get_stock_image_bytes(const std::string& filepath,
                                 ImageFormat format, const RenderOptions& opts)
{
    QFileInfo file_info(filepath.c_str());
    return get_image_bytes<QFileInfo>(file_info, format, opts);
}

namespace
{

template <typename T>
QByteArray get_LiveDesign_image_bytes(const T& input, ImageFormat format,
                                      const RenderOptions& opts)
{
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);
    SketcherModel sketcher_model;
    Scene scene(&mol_model, &sketcher_model);
    init_molviewer_image(mol_model, sketcher_model, input, opts);

    return get_image_bytes<Scene>(scene, format, opts);
}
} // namespace

QByteArray get_LiveDesign_image_bytes(const RDKit::ROMol& mol,
                                      ImageFormat format,
                                      const RenderOptions& opts)
{
    return get_LiveDesign_image_bytes<RDKit::ROMol>(mol, format, opts);
}

QByteArray get_LiveDesign_image_bytes(const RDKit::ChemicalReaction& rxn,
                                      ImageFormat format,
                                      const RenderOptions& opts)
{
    return get_LiveDesign_image_bytes<RDKit::ChemicalReaction>(rxn, format,
                                                               opts);
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
    QUndoStack undo_stack;
    MolModel mol_model(&undo_stack);
    SketcherModel sketcher_model;
    Scene scene(&mol_model, &sketcher_model);
    for (auto rdmol : all_rdmols) {
        init_molviewer_image(mol_model, sketcher_model, *rdmol, opts);
        qreal cur_scale =
            get_scale(scene.getSceneItemsBoundingRect(), opts.width_height);
        best_scale = std::min(best_scale, cur_scale);
        mol_model.clear();
    }
    return best_scale;
}

// SIP has issues wrapping these functions if they're defined using a template,
// so these definitions just call the templated function that contains the
// actual implementation
qreal get_autoscale_value(const RDKit::ROMol& rdmol, const RenderOptions& opts)
{
    return get_image_scale_for_mol_or_rxn(rdmol, opts);
}

qreal get_autoscale_value(const RDKit::ChemicalReaction& rxn,
                          const RenderOptions& opts)
{
    return get_image_scale_for_mol_or_rxn(rxn, opts);
}

} // namespace sketcher
} // namespace schrodinger
