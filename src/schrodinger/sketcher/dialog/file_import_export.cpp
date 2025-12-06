#include "schrodinger/sketcher/dialog/file_import_export.h"

#include <fmt/format.h>
#include <fmt/ranges.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_stream.h"
#include "schrodinger/sketcher/image_generation.h"

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
        {Format::SMILES, "SMILES"},
        {Format::EXTENDED_SMILES, "Extended SMILES"},
        {Format::SMARTS, "SMARTS"},
        {Format::EXTENDED_SMARTS, "Extended SMARTS"},
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

FormatList<Format> get_import_formats()
{

    auto mol_import_formats = get_mol_import_formats();
    auto rxn_import_formats = get_reaction_import_formats();

    FormatList<Format> import_formats;
    for (const auto& [format, label] : mol_import_formats) {
        auto extensions = rdkit_extensions::get_mol_extensions(format);
        // we skip SMARTS and EXTENDED_SMARTS since those don't have any
        // associated extensions
        if (!extensions.empty()) {
            import_formats.push_back({format, label, extensions});
        }
    }
    for (const auto& [format, label] : rxn_import_formats) {
        auto extensions = rdkit_extensions::get_rxn_extensions(format);
        import_formats.push_back({format, label, extensions});
    }
    return import_formats;
}

FormatList<Format> get_standard_export_formats()
{
    std::vector<std::tuple<Format, std::string>> mol_export_formats = {
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
    };

    FormatList<Format> export_formats;
    for (const auto& [format, label] : mol_export_formats) {
        auto extensions = rdkit_extensions::get_mol_extensions(format);
        export_formats.push_back({format, label, extensions});
    }
    return export_formats;
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

    FormatList<Format> export_formats;
    for (const auto& [format, label] : rxn_export_formats) {
        auto extensions = rdkit_extensions::get_rxn_extensions(format);
        export_formats.push_back({format, label, extensions});
    }
    return export_formats;
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
