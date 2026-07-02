/* -------------------------------------------------------------------------
 * Schrodinger Sketcher Application
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/val.h>
#else
#include "crash_handler.h"
#endif

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <QAbstractButton>
#include <QApplication>
#include <QByteArray>
#include <QColor>
#include <QFile>
#include <QHash>
#include <QIcon>
#include <QJsonDocument>
#include <QJsonObject>
#include <QString>
#include <QStyleHints>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/public_constants.h"
#include "schrodinger/sketcher/sketcher_widget.h"

using schrodinger::rdkit_extensions::Format;
using schrodinger::sketcher::CarbonLabels;
using schrodinger::sketcher::ColorScheme;
using schrodinger::sketcher::ImageFormat;
using schrodinger::sketcher::RenderOptions;
using schrodinger::sketcher::SketcherWidget;
using schrodinger::sketcher::StereoLabels;

// For the WebAssembly build, we need to be able to get the sketcher
// instance we are running from a function/static method. We'll use a
// singleton for that.
SketcherWidget& get_sketcher_instance()
{
    static SketcherWidget instance;
    return instance;
}

void sketcher_import_text(const std::string& text)
{
    auto& sk = get_sketcher_instance();
    sk.addFromString(text);
}

std::string sketcher_export_text(Format format)
{
    auto& sk = get_sketcher_instance();
    return sk.getString(format);
}

std::string sketcher_export_image(ImageFormat format)
{
    auto& sk = get_sketcher_instance();
    return sk.getImageBytes(format).toBase64().toStdString();
}

#ifdef __EMSCRIPTEN__
// code required for the get_image_bytes wrapping
namespace
{

emscripten::val qbyte_array_to_uint8_array(QByteArray bytes)
{
    auto byte_view = emscripten::val(emscripten::typed_memory_view(
        static_cast<std::size_t>(bytes.size()),
        reinterpret_cast<unsigned char*>(bytes.data())));
    return emscripten::val::global("Uint8Array").new_(byte_view);
}

bool has_option(const emscripten::val& options, const char* key)
{
    return options.hasOwnProperty(key);
}

int parse_js_index(const std::string& index, const char* option_name)
{
    std::size_t parsed_chars = 0;
    int parsed_index = 0;
    try {
        parsed_index = std::stoi(index, &parsed_chars);
    } catch (const std::exception&) {
        throw std::invalid_argument("RenderOptions." +
                                    std::string(option_name) +
                                    " contains invalid index '" + index + "'");
    }
    if (parsed_chars != index.size() || parsed_index < 0) {
        throw std::invalid_argument("RenderOptions." +
                                    std::string(option_name) +
                                    " contains invalid index '" + index + "'");
    }
    return parsed_index;
}

QColor parse_js_color(const emscripten::val& color_val, const char* option_name)
{
    auto color_string = color_val.as<std::string>();
    if (color_string == "transparent") {
        return Qt::transparent;
    }

    auto color = QColor::fromString(QString::fromStdString(color_string));
    if (!color.isValid()) {
        throw std::invalid_argument(
            "RenderOptions." + std::string(option_name) +
            " contains invalid color '" + color_string + "'");
    }
    return color;
}

std::vector<std::string> get_js_object_keys(const emscripten::val& object)
{
    auto keys =
        emscripten::val::global("Object").call<emscripten::val>("keys", object);
    return emscripten::vecFromJSArray<std::string>(keys);
}

void assert_js_object(const emscripten::val& object, const char* option_name)
{
    if (object.typeOf().as<std::string>() != "object" || !object.as<bool>()) {
        throw std::invalid_argument(
            "RenderOptions." + std::string(option_name) + " must be an object");
    }
}

void read_string_hash_option(const emscripten::val& options,
                             const char* option_name,
                             QHash<int, std::string>& target)
{
    if (!has_option(options, option_name)) {
        return;
    }

    auto js_hash = options[option_name];
    assert_js_object(js_hash, option_name);
    for (const auto& key : get_js_object_keys(js_hash)) {
        target[parse_js_index(key, option_name)] =
            js_hash[key].as<std::string>();
    }
}

void read_color_hash_option(const emscripten::val& options,
                            const char* option_name, QHash<int, QColor>& target)
{
    if (!has_option(options, option_name)) {
        return;
    }

    auto js_hash = options[option_name];
    assert_js_object(js_hash, option_name);
    for (const auto& key : get_js_object_keys(js_hash)) {
        target[parse_js_index(key, option_name)] =
            parse_js_color(js_hash[key], option_name);
    }
}

void read_size_option(const emscripten::val& options, RenderOptions& opts)
{
    if (has_option(options, "width_height")) {
        auto js_size = options["width_height"];
        assert_js_object(js_size, "width_height");
        if (has_option(js_size, "width") && has_option(js_size, "height")) {
            opts.width_height =
                QSize(js_size["width"].as<int>(), js_size["height"].as<int>());
        } else if (has_option(js_size, "0") && has_option(js_size, "1")) {
            opts.width_height =
                QSize(js_size[0].as<int>(), js_size[1].as<int>());
        } else {
            throw std::invalid_argument(
                "RenderOptions.width_height must provide width and height");
        }
    }
    if (has_option(options, "width")) {
        opts.width_height.setWidth(options["width"].as<int>());
    }
    if (has_option(options, "height")) {
        opts.width_height.setHeight(options["height"].as<int>());
    }
}

RenderOptions render_options_from_js(const emscripten::val& options)
{
    RenderOptions opts;
    if (options.typeOf().as<std::string>() == "undefined") {
        return opts;
    }
    assert_js_object(options, "options");

    read_size_option(options, opts);
    if (has_option(options, "background_color")) {
        opts.background_color =
            parse_js_color(options["background_color"], "background_color");
    }
    if (has_option(options, "scale")) {
        opts.scale = options["scale"].as<qreal>();
    }
    if (has_option(options, "trim_image")) {
        opts.trim_image = options["trim_image"].as<bool>();
    }
    if (has_option(options, "font_size")) {
        opts.font_size = options["font_size"].as<int>();
    }
    if (has_option(options, "bond_width_scale")) {
        opts.bond_width_scale = options["bond_width_scale"].as<qreal>();
    }
    read_string_hash_option(options, "rdatom_index_to_label",
                            opts.rdatom_index_to_label);
    read_color_hash_option(options, "rdatom_index_to_halo_color",
                           opts.rdatom_index_to_halo_color);
    read_color_hash_option(options, "rdbond_index_to_halo_color",
                           opts.rdbond_index_to_halo_color);
    read_color_hash_option(options, "rdatom_index_to_line_color",
                           opts.rdatom_index_to_line_color);
    read_color_hash_option(options, "rdbond_index_to_line_color",
                           opts.rdbond_index_to_line_color);
    if (has_option(options, "show_stereo_annotations")) {
        opts.show_stereo_annotations =
            options["show_stereo_annotations"].as<StereoLabels>();
    }
    if (has_option(options, "show_absolute_stereo_groups")) {
        opts.show_absolute_stereo_groups =
            options["show_absolute_stereo_groups"].as<bool>();
    }
    if (has_option(options, "show_simplified_stereo_annotation")) {
        opts.show_simplified_stereo_annotation =
            options["show_simplified_stereo_annotation"].as<bool>();
    }
    if (has_option(options, "show_symbol_for_H_isotopes")) {
        opts.show_symbol_for_H_isotopes =
            options["show_symbol_for_H_isotopes"].as<bool>();
    }
    if (has_option(options, "carbon_labels")) {
        opts.carbon_labels = options["carbon_labels"].as<CarbonLabels>();
    }
    if (has_option(options, "color_scheme")) {
        opts.color_scheme = options["color_scheme"].as<ColorScheme>();
    }
    return opts;
}

} // namespace

emscripten::val get_image_bytes_from_text(const std::string& text,
                                          ImageFormat format)
{
    return qbyte_array_to_uint8_array(
        schrodinger::sketcher::get_image_bytes(text, format));
}

emscripten::val get_image_bytes_from_text(const std::string& text,
                                          ImageFormat format,
                                          const emscripten::val& options)
{
    return qbyte_array_to_uint8_array(schrodinger::sketcher::get_image_bytes(
        text, format, render_options_from_js(options)));
}
#endif

void sketcher_clear()
{
    auto& sk = get_sketcher_instance();
    sk.clear();
}

bool sketcher_is_empty()
{
    auto& sk = get_sketcher_instance();
    return sk.isEmpty();
}

bool sketcher_has_monomers()
{
    auto& sk = get_sketcher_instance();
    auto mol = sk.getRDKitMolecule();
    return schrodinger::rdkit_extensions::isMonomeric(*mol);
}

// Retained as a no-op for backwards compatibility with external callers;
// ATOMISTIC_OR_MONOMERIC is now the default interface type (SKETCH-2735).
void sketcher_allow_monomeric(bool /* allow_monomeric */)
{
}

void sketcher_load_custom_monomers(const std::string& json)
{
    auto& db = schrodinger::rdkit_extensions::MonomerDatabase::instance();
    db.loadMonomersFromJson(json);
}

void sketcher_load_custom_monomers_from_sql(const std::string& sql)
{
    auto& db = schrodinger::rdkit_extensions::MonomerDatabase::instance();
    db.loadMonomersFromSql(sql);
}

void sketcher_insert_custom_monomers(const std::string& json)
{
    auto& db = schrodinger::rdkit_extensions::MonomerDatabase::instance();
    db.insertMonomersFromJson(json);
}

void sketcher_reset_custom_monomers()
{
    auto& db = schrodinger::rdkit_extensions::MonomerDatabase::instance();
    db.resetMonomerDefinitions();
}

/**
 * Programmatically click a button by its Qt objectName.
 * Used by e2e tests to interact with popup buttons that can't be targeted
 * via browser events in the WASM environment.
 *
 * Throws std::runtime_error if no button with the given name is found,
 * which typically indicates the test needs to be updated.
 */
void sketcher_click_button(const std::string& name)
{
    auto& sk = get_sketcher_instance();
    auto* button = sk.findChild<QAbstractButton*>(QString::fromStdString(name));
    if (!button) {
        throw std::runtime_error(
            "sketcher_click_button: no button found with objectName '" + name +
            "'");
    }
    button->click();
}

/**
 * Return the position and size of any child widget found by its Qt
 * objectName. Returns a JSON object with {x, y, width, height}.
 * Coordinates are relative to the sketcher widget's top-left corner.
 * Returns "{}" if no widget with the given name is found.
 */
std::string sketcher_get_widget_rect(const std::string& object_name)
{
    auto& sk = get_sketcher_instance();
    const QString name = QString::fromStdString(object_name);
    QWidget* widget = nullptr;
    const auto candidates = sk.findChildren<QWidget*>(name);
    for (auto* w : candidates) {
        if (w->isVisible()) {
            widget = w;
            break;
        }
    }
    if (!widget) {
        return "{}";
    }
    const QPoint topLeft = widget->mapTo(&sk, QPoint(0, 0));
    QJsonObject result;
    result["x"] = topLeft.x();
    result["y"] = topLeft.y();
    result["width"] = widget->width();
    result["height"] = widget->height();
    return QJsonDocument(result).toJson(QJsonDocument::Compact).toStdString();
}

void sketcher_changed()
{
#ifdef __EMSCRIPTEN__
    EM_ASM({
        if (Module.sketcher_changed_callback) {
            setTimeout(
                function() {
                    if (Module.sketcher_changed_callback) {
                        Module.sketcher_changed_callback();
                    }
                },
                100);
        }
    });
#endif
}

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(sketcher)
{
    emscripten::enum_<Format>("Format")
        .value("AUTO_DETECT", Format::AUTO_DETECT)
        .value("RDMOL_BINARY_BASE64", Format::RDMOL_BINARY_BASE64)
        .value("SMILES", Format::SMILES)
        .value("EXTENDED_SMILES", Format::EXTENDED_SMILES)
        .value("SMARTS", Format::SMARTS)
        .value("EXTENDED_SMARTS", Format::EXTENDED_SMARTS)
        .value("MDL_MOLV2000", Format::MDL_MOLV2000)
        .value("MDL_MOLV3000", Format::MDL_MOLV3000)
        .value("MAESTRO", Format::MAESTRO)
        .value("INCHI", Format::INCHI)
        .value("INCHI_KEY", Format::INCHI_KEY)
        .value("PDB", Format::PDB)
        .value("MOL2", Format::MOL2)
        .value("XYZ", Format::XYZ)
        .value("MRV", Format::MRV)
        .value("CDXML", Format::CDXML)
        .value("HELM", Format::HELM)
        .value("FASTA_PEPTIDE", Format::FASTA_PEPTIDE)
        .value("FASTA_DNA", Format::FASTA_DNA)
        .value("FASTA_RNA", Format::FASTA_RNA)
        .value("FASTA", Format::FASTA);

    emscripten::enum_<ImageFormat>("ImageFormat")
        .value("PNG", ImageFormat::PNG)
        .value("SVG", ImageFormat::SVG);

    emscripten::enum_<StereoLabels>("StereoLabels")
        .value("NONE", StereoLabels::NONE)
        .value("KNOWN", StereoLabels::KNOWN)
        .value("ALL", StereoLabels::ALL);

    emscripten::enum_<CarbonLabels>("CarbonLabels")
        .value("NONE", CarbonLabels::NONE)
        .value("TERMINAL", CarbonLabels::TERMINAL)
        .value("ALL", CarbonLabels::ALL);

    emscripten::enum_<ColorScheme>("ColorScheme")
        .value("DEFAULT", ColorScheme::DEFAULT)
        .value("AVALON", ColorScheme::AVALON)
        .value("CDK", ColorScheme::CDK)
        .value("DARK_MODE", ColorScheme::DARK_MODE)
        .value("BLACK_WHITE", ColorScheme::BLACK_WHITE)
        .value("WHITE_BLACK", ColorScheme::WHITE_BLACK);

    emscripten::function("sketcher_import_text", &sketcher_import_text);
    emscripten::function("sketcher_export_text", &sketcher_export_text);
    emscripten::function("sketcher_export_image", &sketcher_export_image);
    emscripten::function(
        "get_image_bytes",
        emscripten::select_overload<emscripten::val(
            const std::string&, ImageFormat)>(&get_image_bytes_from_text));
    emscripten::function(
        "get_image_bytes",
        emscripten::select_overload<emscripten::val(
            const std::string&, ImageFormat, const emscripten::val&)>(
            &get_image_bytes_from_text));
    emscripten::function("sketcher_clear", &sketcher_clear);
    emscripten::function("sketcher_is_empty", &sketcher_is_empty);
    emscripten::function("sketcher_has_monomers", &sketcher_has_monomers);
    emscripten::function("sketcher_allow_monomeric", &sketcher_allow_monomeric);
    emscripten::function("sketcher_load_custom_monomers",
                         &sketcher_load_custom_monomers);
    emscripten::function("sketcher_load_custom_monomers_from_sql",
                         &sketcher_load_custom_monomers_from_sql);
    emscripten::function("sketcher_insert_custom_monomers",
                         &sketcher_insert_custom_monomers);
    emscripten::function("sketcher_reset_custom_monomers",
                         &sketcher_reset_custom_monomers);
    emscripten::function("_sketcher_get_widget_rect",
                         &sketcher_get_widget_rect);
    emscripten::function("_sketcher_click_button", &sketcher_click_button);
    // see sketcher_changed_callback above
}
#endif

void apply_stylesheet(QApplication& app)
{
    // In Qt 6.8 and newer, Qt will try to automatically apply a dark mode color
    // scheme if the system and/or browser is set to dark mode. The result looks
    // terrible, so switch back to light mode.
#if QT_VERSION >= QT_VERSION_CHECK(6, 8, 0)
    QApplication::styleHints()->setColorScheme(Qt::ColorScheme::Light);
#endif

    QFile styleFile(":resources/schrodinger_livedesign.qss");
    bool success = styleFile.open(QFile::ReadOnly);
    if (!success) {
        throw std::runtime_error("Could not open style sheet file");
    }
    QString style(styleFile.readAll());
    app.setStyleSheet(style);
}

int main(int argc, char** argv)
{
#ifndef __EMSCRIPTEN__
    schrodinger::install_crash_handlers();
#endif

    QApplication application(argc, argv);
#ifdef SKETCHER_STATIC_DEFINE
    Q_INIT_RESOURCE(sketcher);
#endif
    QApplication::setWindowIcon(QIcon(":icons/sketcher-logo.svg"));

#ifdef __EMSCRIPTEN__
    // Only apply this stylesheet for the WASM build
    apply_stylesheet(application);
    auto& sk = get_sketcher_instance();
    // Qt::WA_AlwaysShowToolTips works around QTBUG-94583
    // (https://bugreports.qt.io/browse/QTBUG-94583), which would otherwise
    // prevent tooltips from showing up. (See SKETCH-2565.)
    sk.setAttribute(Qt::WA_AlwaysShowToolTips);
    QObject::connect(&sk, &SketcherWidget::moleculeChanged, &sketcher_changed);
    QObject::connect(&sk, &SketcherWidget::representationChanged,
                     &sketcher_changed);
#else
    SketcherWidget sk;
#endif

    sk.show();
    return application.exec();
}
