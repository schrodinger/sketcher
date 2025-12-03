#pragma once

#include <tuple>

#include <QString>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{

namespace sketcher
{

enum class ImageFormat;

// Collection specifying the relationship for permitted formats via tuples
// containing (format enum, menu label, allowable extensions). These are
// stored as an ordered list to preserve menu initialization order.
template <class T> using FormatList =
    std::vector<std::tuple<T, std::string, std::vector<std::string>>>;

/**
 * @return a list of all atomistic import formats with a human readable name for
 * each
 */
SKETCHER_API std::vector<std::tuple<rdkit_extensions::Format, std::string>>
get_mol_import_formats();

/**
 * @return a list of all monomeric import formats with a human readable name for
 * each
 */
SKETCHER_API std::vector<std::tuple<rdkit_extensions::Format, std::string>>
get_monomoric_import_formats();

/**
 * @return a list of all atomistic reaction import formats with a human readable
 * name for each
 */
SKETCHER_API std::vector<std::tuple<rdkit_extensions::Format, std::string>>
get_reaction_import_formats();

/**
 * @return list of importable (format enum, menu label, allowable extensions)
 */
SKETCHER_API FormatList<rdkit_extensions::Format> get_import_formats();

/**
 * @return list of exportable (format enum, menu label, allowable extensions)
 *  for both standard molecules, reactions, and image formats
 */
SKETCHER_API FormatList<rdkit_extensions::Format> get_standard_export_formats();
SKETCHER_API FormatList<rdkit_extensions::Format> get_reaction_export_formats();
SKETCHER_API FormatList<ImageFormat> get_image_export_formats();

/**
 * @param file_path file to read
 * @return full contents of that file as a single string
 */
SKETCHER_API std::string get_file_text(const std::string& file_path);

/**
 * @param label filter label prefix
 * @param extensions filter extensions to append
 * @return name filter to use for the import/export file dialogs
 */
SKETCHER_API QString get_filter_name(
    const std::string& label, const std::vector<std::string>& extensions);

} // namespace sketcher
} // namespace schrodinger
