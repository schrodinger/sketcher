#ifdef __EMSCRIPTEN__

#include "image_generation_from_js.h"

#include <emscripten/bind.h>

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include <QColor>
#include <QHash>
#include <QString>
#include <QSize>

using schrodinger::sketcher::CarbonLabels;
using schrodinger::sketcher::ColorScheme;
using schrodinger::sketcher::RenderOptions;
using schrodinger::sketcher::StereoLabels;

static bool has_option(const emscripten::val& options, const char* key)
{
    return options.hasOwnProperty(key);
}

static int parse_js_index(const std::string& index, const char* option_name)
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

static QColor parse_js_color(const emscripten::val& color_val,
                             const char* option_name)
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

static std::vector<std::string>
get_js_object_keys(const emscripten::val& object)
{
    auto keys =
        emscripten::val::global("Object").call<emscripten::val>("keys", object);
    return emscripten::vecFromJSArray<std::string>(keys);
}

static void assert_js_object(const emscripten::val& object,
                             const char* option_name)
{
    if (object.typeOf().as<std::string>() != "object" || !object.as<bool>()) {
        throw std::invalid_argument(
            "RenderOptions." + std::string(option_name) + " must be an object");
    }
}

static void read_string_hash_option(const emscripten::val& options,
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

static void read_color_hash_option(const emscripten::val& options,
                                   const char* option_name,
                                   QHash<int, QColor>& target)
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

static void read_size_option(const emscripten::val& options,
                             RenderOptions& opts)
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

#endif
