#pragma once

#include <tuple>

#include <QString>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/public_constants.h"

namespace schrodinger
{

namespace sketcher
{

enum class ImageFormat;
enum class MoleculeType;
enum class ToolSet;

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
 * Determine which kinds of structure the user is allowed to import right now.
 *
 * @param interface_type which kinds of structure the interface supports
 * @param cur_mol_type what the Sketcher currently contains
 * @param replace_content whether an import replaces the current contents
 * @return the permitted molecule types, as an InterfaceType bitmask
 */
SKETCHER_API InterfaceTypeType get_importable_mol_types(
    const InterfaceTypeType interface_type, const MoleculeType cur_mol_type,
    const bool replace_content);

/**
 * @param interface_type which kinds of structure the interface supports
 * @param cur_mol_type what the Sketcher currently contains
 * @param replace_content whether an import replaces the current contents
 * @param tool_set which tab the user is on
 * @return list of importable (format enum, menu label, allowable extensions).
 * Sequence formats are listed first on the Monomer tab and last on the
 * Atomistic tab, and are omitted entirely when they can't be imported.
 */
SKETCHER_API FormatList<rdkit_extensions::Format>
get_import_formats(const InterfaceTypeType interface_type,
                   const MoleculeType cur_mol_type, const bool replace_content,
                   const ToolSet tool_set);

/**
 * File extensions can't distinguish the FASTA sub-formats, so
 * rdkit_extensions::get_file_format() reports the unreadable Format::FASTA for
 * all of them.
 *
 * @param format the format reported for a file
 * @return a format that can actually be read: unchanged unless the input was
 * Format::FASTA, in which case Format::FASTA_PEPTIDE
 */
SKETCHER_API rdkit_extensions::Format
resolve_ambiguous_import_format(const rdkit_extensions::Format format);

/**
 * @return the molecule, sequence, and reaction export formats with a human
 * readable name for each. Callers that write to a file should prefer
 * get_standard_export_formats() / get_reaction_export_formats(), which also
 * report extensions; these are for destinations with no filename, such as the
 * clipboard, where compression doesn't apply.
 */
SKETCHER_API std::vector<std::tuple<rdkit_extensions::Format, std::string>>
get_mol_and_seq_export_formats();
SKETCHER_API std::vector<std::tuple<rdkit_extensions::Format, std::string>>
get_rxn_export_formats();

/**
 * @return list of exportable (format enum, menu label, allowable extensions)
 *  for both standard molecules, reactions, and image formats. Each format's
 *  compressed extensions get their own "[compressed]" entry directly below the
 *  uncompressed one.
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
