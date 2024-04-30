#pragma once

#include <tuple>

#include <QFileInfo>
#include <QList>
#include <QString>
#include <QStringList>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace rdkit_extensions
{
enum class Format;
} // namespace rdkit_extensions

namespace sketcher
{
enum class ImageFormat;

// Collection specifying the relationship for permitted formats via tuples
// containing (format enum, menu label, allowable extensions). These are
// stored as an ordered list to preserve menu initialization order.
template <class T> using FormatList =
    QList<std::tuple<T, QString, QStringList>>;

/**
 * All supported export labels for both standard molecules and reactions
 */
extern SKETCHER_API const FormatList<rdkit_extensions::Format> STANDARD_FORMATS;
extern SKETCHER_API const FormatList<rdkit_extensions::Format> REACTION_FORMATS;
extern SKETCHER_API const FormatList<ImageFormat> IMAGE_FORMATS;

/**
 * @param file_path file to read
 * @return full contents of that file as a single string
 */
SKETCHER_API std::string get_file_text(const std::string& file_path);

/**
 * @return concatenated name filters to use for the import file dialog
 */
SKETCHER_API QString get_import_name_filters();

/**
 * @param data [format, label, extensions] to iterate over
 * @return name filters to use for the import/export file dialogs
 */
template <class T>
QList<std::tuple<T, QString>> get_name_filters(const FormatList<T>& format_list)
{
    QList<std::tuple<T, QString>> filters;
    for (const auto& [format, label, extensions] : format_list) {
        if (!extensions.isEmpty()) {
            auto filter = label + " (*" + extensions.join(" *") + ")";
            filters.append({format, filter});
        }
    }
    return filters;
}

/**
 * @param data [format, label, extensions] to iterate over
 * @param request_format the requested format
 * @return supported file extensions for the given format
 */
template <class T> QStringList
get_file_extensions(const FormatList<T>& format_list, const T& request_format)
{
    for (const auto& [format, _, extensions] : format_list) {
        if (format == request_format) {
            return extensions;
        }
    }
    throw std::runtime_error("Unknown format requested: " +
                             std::to_string(static_cast<int>(request_format)));
}

} // namespace sketcher
} // namespace schrodinger
