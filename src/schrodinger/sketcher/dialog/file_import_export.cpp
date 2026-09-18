#include "schrodinger/sketcher/dialog/file_import_export.h"

#include <boost/algorithm/string/predicate.hpp>
#include <fmt/format.h>
#include <fmt/ranges.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_stream.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

std::vector<std::tuple<Format, std::string>> get_mol_import_formats()
{
    return {
        {Format::MDL_MOLV3000, "MDL SD"},
        {Format::MAESTRO, "Maestro"},
        // On input there is no difference between SMILES and CXSMILES (or
        // between SMARTS and CXSMARTS); both are handled by the same parser,
        // so we offer a single entry for each rather than one per variant.
        {Format::SMILES, "SMILES"},
        {Format::SMARTS, "SMARTS"},
        {Format::INCHI, "InChI"},
        {Format::MOL2, "MOL2"},
        {Format::PDB, "PDB"},
        {Format::XYZ, "XYZ"},
        {Format::MRV, "Marvin Document"},
        {Format::CDXML, "ChemDraw XML"},
    };
}

std::vector<std::tuple<Format, std::string>> get_reaction_import_formats()
{
    return {
        {Format::MDL_MOLV3000, "MDL RXN"},
        {Format::SMILES, "Reaction SMILES"},
    };
}

std::vector<std::tuple<Format, std::string>> get_monomoric_import_formats()
{
    return {
        {Format::HELM, "HELM"},
        {Format::FASTA_PEPTIDE, "FASTA Peptide"},
        {Format::FASTA_DNA, "FASTA DNA"},
        {Format::FASTA_RNA, "FASTA RNA"},
    };
}

namespace
{

/**
 * @return the sequence formats to offer in a file dialog, with a human
 * readable name for each. This differs from get_monomoric_import_formats(),
 * which lists the FASTA sub-formats used when the caller knows what kind of
 * sequence it has; extensions can't distinguish them, so a file dialog offers
 * a single FASTA entry and resolve_ambiguous_import_format() picks a reader.
 */
std::vector<std::tuple<Format, std::string>> get_seq_file_import_formats()
{
    return {
        {Format::HELM, "HELM"},
        {Format::FASTA, "FASTA"},
    };
}

/**
 * @param ext a file extension, including the leading dot
 * @return whether the extension denotes a compressed file
 */
bool is_compressed_extension(const std::string& ext)
{
    return boost::algorithm::iends_with(ext, "gz") ||
           boost::algorithm::iends_with(ext, "zst");
}

/**
 * @param format the format to look up extensions for
 * @return the molecule extensions for that format. The SMILES entry stands in
 * for CXSMILES as well (the two are read and written by the same parser), so
 * it also claims the CXSMILES extensions.
 */
std::vector<std::string> get_mol_import_extensions(const Format format)
{
    auto extensions = rdkit_extensions::get_mol_extensions(format);
    if (format == Format::SMILES) {
        auto cxsmiles_extensions =
            rdkit_extensions::get_mol_extensions(Format::EXTENDED_SMILES);
        extensions.insert(extensions.end(), cxsmiles_extensions.begin(),
                          cxsmiles_extensions.end());
    }
    return extensions;
}

/**
 * Build menu entries for the given formats, splitting each format's extensions
 * into an uncompressed entry and a "[compressed]" entry so that no single row
 * lists both. Formats with no extensions at all are omitted, since they can't
 * be selected in a file dialog.
 *
 * @param formats (format enum, menu label) pairs, in menu order
 * @param get_extensions callable returning the extensions for one format
 * @param[out] uncompressed_entries uncompressed entries are appended here
 * @param[out] compressed_entries compressed entries are appended here
 */
template <class T, class F> void
append_split_entries(const std::vector<std::tuple<T, std::string>>& formats,
                     F get_extensions, FormatList<T>& uncompressed_entries,
                     FormatList<T>& compressed_entries)
{
    for (const auto& [format, label] : formats) {
        std::vector<std::string> uncompressed;
        std::vector<std::string> compressed;
        for (const auto& ext : get_extensions(format)) {
            (is_compressed_extension(ext) ? compressed : uncompressed)
                .push_back(ext);
        }
        if (!uncompressed.empty()) {
            uncompressed_entries.push_back({format, label, uncompressed});
        }
        if (!compressed.empty()) {
            compressed_entries.push_back(
                {format, label + " [compressed]", compressed});
        }
    }
}

/**
 * @return the given formats as menu entries, with all uncompressed entries
 * first and all "[compressed]" entries after them
 */
template <class T, class F> FormatList<T>
split_compressed_formats(const std::vector<std::tuple<T, std::string>>& formats,
                         F get_extensions)
{
    FormatList<T> entries;
    FormatList<T> compressed_entries;
    append_split_entries(formats, get_extensions, entries, compressed_entries);
    entries.insert(entries.end(), compressed_entries.begin(),
                   compressed_entries.end());
    return entries;
}

} // namespace

InterfaceTypeType
get_importable_mol_types(const InterfaceTypeType interface_type,
                         const MoleculeType cur_mol_type,
                         const bool replace_content)
{
    // If the interface only offers one kind of structure, that's all we can
    // import. Otherwise, unless the import is going to replace what's already
    // in the Sketcher, we're limited to whatever is already there.
    if (interface_type == InterfaceType::ATOMISTIC ||
        (!replace_content && cur_mol_type == MoleculeType::ATOMISTIC)) {
        return InterfaceType::ATOMISTIC;
    }
    if (interface_type == InterfaceType::MONOMERIC ||
        (!replace_content && cur_mol_type == MoleculeType::MONOMERIC)) {
        return InterfaceType::MONOMERIC;
    }
    return InterfaceType::ATOMISTIC_OR_MONOMERIC;
}

FormatList<Format> get_import_formats(const InterfaceTypeType interface_type,
                                      const MoleculeType cur_mol_type,
                                      const bool replace_content)
{
    auto allowed_mol_types =
        get_importable_mol_types(interface_type, cur_mol_type, replace_content);

    FormatList<Format> atomistic_entries;
    FormatList<Format> atomistic_compressed;
    if (allowed_mol_types & InterfaceType::ATOMISTIC) {
        append_split_entries(get_mol_import_formats(),
                             get_mol_import_extensions, atomistic_entries,
                             atomistic_compressed);
        append_split_entries(get_reaction_import_formats(),
                             rdkit_extensions::get_rxn_extensions,
                             atomistic_entries, atomistic_compressed);
        atomistic_entries.insert(atomistic_entries.end(),
                                 atomistic_compressed.begin(),
                                 atomistic_compressed.end());
    }

    FormatList<Format> monomeric_entries;
    if (allowed_mol_types & InterfaceType::MONOMERIC) {
        monomeric_entries =
            split_compressed_formats(get_seq_file_import_formats(),
                                     rdkit_extensions::get_seq_extensions);
    }

    // Laura Beck (SKETCH-2516): on the Monomer tab the sequence formats belong
    // at the top, since everything else is atom-based; elsewhere they go at the
    // bottom. A visual separator between the two groups isn't possible here,
    // because QFileDialog::getOpenFileContent only accepts a filter string.
    FormatList<Format> import_formats;
    if (interface_type == InterfaceType::MONOMERIC) {
        import_formats = std::move(monomeric_entries);
        import_formats.insert(import_formats.end(), atomistic_entries.begin(),
                              atomistic_entries.end());
    } else {
        import_formats = std::move(atomistic_entries);
        import_formats.insert(import_formats.end(), monomeric_entries.begin(),
                              monomeric_entries.end());
    }
    return import_formats;
}

FormatList<Format> get_standard_export_formats()
{
    std::vector<std::tuple<Format, std::string>> mol_and_seq_export_formats = {
        // Forbid MDL_MOLV2000 on export; potential stereo ambiguities
        {Format::MDL_MOLV3000, "MDL SD V3000"},
        {Format::MAESTRO, "Maestro"},
        {Format::SMILES, "SMILES"},
        {Format::EXTENDED_SMILES, "Extended SMILES"},
        {Format::SMARTS, "SMARTS"},
        {Format::EXTENDED_SMARTS, "Extended SMARTS"},
        {Format::INCHI, "InChI"},
        {Format::INCHI_KEY, "InChIKey"},
        {Format::PDB, "PDB"},
        {Format::XYZ, "XYZ"},
        {Format::MRV, "Marvin Document"},
        // Sequence formats are always offered; rdkit_extensions::to_string
        // handles atomistic <-> monomeric conversion on the fly, and conversion
        // failures surface through the dialog's existing error path.
        {Format::HELM, "HELM"},
        {Format::FASTA, "FASTA"},
    };

    return split_compressed_formats(
        mol_and_seq_export_formats,
        rdkit_extensions::get_mol_and_seq_extensions);
};

FormatList<Format> get_reaction_export_formats()
{
    std::vector<std::tuple<Format, std::string>> rxn_export_formats = {
        // Forbid MDL_MOLV2000 on export; potential stereo ambiguities
        {Format::MDL_MOLV3000, "MDL RXN V3000"},
        {Format::SMILES, "Reaction SMILES"},
        {Format::EXTENDED_SMILES, "Extended Reaction SMILES"},
        {Format::SMARTS, "Reaction SMARTS"},
        {Format::EXTENDED_SMARTS, "Extended Reaction SMARTS"},
    };

    return split_compressed_formats(rxn_export_formats,
                                    rdkit_extensions::get_rxn_extensions);
};

// We define get_image_formats as a function for consistency with
// the above, even though it doesn't depend on any other values
FormatList<ImageFormat> get_image_export_formats()
{
    std::vector<std::tuple<ImageFormat, std::string>> image_export_formats = {
        {ImageFormat::PNG, "PNG"},
        {ImageFormat::SVG, "SVG"},
    };

    FormatList<ImageFormat> export_formats;
    for (const auto& [format, label] : image_export_formats) {
        auto extension = get_image_extension(format);
        export_formats.push_back({format, label, {extension}});
    }
    return export_formats;
}

Format resolve_ambiguous_import_format(const Format format)
{
    // Every FASTA extension maps to Format::FASTA, which can't be read
    // directly. The three sub-formats share an alphabet, so there's no way to
    // tell them apart from the contents; assume peptide, as AUTO_DETECT_FORMATS
    // does.
    return format == Format::FASTA ? Format::FASTA_PEPTIDE : format;
}

std::string get_file_text(const std::string& file_path)
{
    rdkit_extensions::maybe_compressed_istream file(file_path);
    if (file.fail()) {
        throw std::runtime_error("Cannot open the file: " + file_path);
    }
    std::string text((std::istreambuf_iterator<char>(file)),
                     std::istreambuf_iterator<char>());
    return text;
}

QString get_filter_name(const std::string& label,
                        const std::vector<std::string>& extensions)
{
    return fmt::format("{} (*{})", label, fmt::join(extensions, " *")).c_str();
}

} // namespace sketcher
} // namespace schrodinger
