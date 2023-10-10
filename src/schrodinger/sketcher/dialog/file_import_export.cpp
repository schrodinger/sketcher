#include "schrodinger/sketcher/dialog/file_import_export.h"

#include <fstream>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/image_generation.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

// SKETCH-1453: Forbid MDL_MOLV2000 on export; potential stereo ambiguities

// File extensions are drawn from Open Babel:
// https://openbabel.org/docs/current/FileFormats/Common_cheminformatics_Formats.html

const FormatList<Format> STANDARD_FORMATS{
    {Format::MDL_MOLV3000, "MDL SD V3000", {".sdf", ".sd", ".mol", ".mdl"}},
    {Format::SMILES, "SMILES", {".smi", ".smiles"}},
    {Format::EXTENDED_SMILES, "Extended SMILES", {".cxsmi", ".cxsmiles"}},
    {Format::SMARTS, "SMARTS", {}},
    {Format::INCHI, "InChI", {".inchi"}},
    {Format::INCHI_KEY, "InChIKey", {}},
    {Format::PDB, "PDB", {".pdb", ".ent"}},
    {Format::XYZ, "XYZ", {".xyz"}},
};

namespace
{
// SHARED-9525: Only allow MAE_FORMAT as an allowable import format for now;
// we exclude it from STANDARD_FORMATS so that it cannot be exported until
// the RDKit MaeWriter can write stereochemistry.
const FormatList<Format> MAE_FORMAT{
    {Format::MAESTRO, "Maestro", {".mae", ".maegz", ".mae.gz"}}};
const FormatList<Format> IMPORT_FORMATS = STANDARD_FORMATS + MAE_FORMAT;
} // unnamed namespace

const FormatList<Format> REACTION_FORMATS{
    {Format::MDL_MOLV3000, "MDL RXN V3000", {".rxn"}},
    {Format::SMILES, "Reaction SMILES", {".rsmi"}},
    {Format::SMARTS, "Reaction SMARTS", {}},
};

const FormatList<ImageFormat> IMAGE_FORMATS{
    {ImageFormat::PNG, "PNG", {".png"}},
    {ImageFormat::SVG, "SVG", {".svg"}},
};

std::string get_file_text(const std::string& file_path)
{
    std::ifstream file(file_path);
    if (file.fail()) {
        throw std::runtime_error("Cannot open the file: " + file_path);
    }
    std::string text((std::istreambuf_iterator<char>(file)),
                     std::istreambuf_iterator<char>());
    return text;
}

Format get_file_format(const QString& file_path)
{
    auto ext = "." + QFileInfo(file_path).completeSuffix();
    for (const auto& [format, _, extensions] : IMPORT_FORMATS) {
        if (extensions.contains(ext)) {
            return format;
        }
    }
    for (const auto& [format, _, extensions] : REACTION_FORMATS) {
        if (extensions.contains(ext)) {
            return format;
        }
    }

    throw std::runtime_error("Unknown file extension: " +
                             file_path.toStdString());
}

QString get_import_name_filters()
{
    QStringList filters;
    for (const auto& format_list : {IMPORT_FORMATS, REACTION_FORMATS}) {
        for (const auto& [_, filter] : get_name_filters(format_list)) {
            filters.append(filter);
        }
    }

    // While we enforce MDL export via V3000 (see SKETCH-1453), we allow import
    // from either v2000 and v3000. Remove any V3000 label to avoid confusion.
    filters.replaceInStrings(" V3000", "");
    return filters.join(";;");
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/dialog/file_import_export.moc"
