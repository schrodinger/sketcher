#include "schrodinger/sketcher/dialog/file_import_export.h"

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_stream.h"
#include "schrodinger/sketcher/image_generation.h"

using ::schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

namespace
{

QStringList to_qstringlist(const std::vector<std::string>& vec)
{
    QStringList list;
    for (const auto& str : vec) {
        list.append(QString::fromStdString(str));
    }
    return list;
}

std::tuple<Format, QString, QStringList> mol_data(Format format,
                                                  const QString& name)
{
    return {format, name,
            to_qstringlist(rdkit_extensions::get_mol_extensions(format))};
}

std::tuple<Format, QString, QStringList> rxn_data(Format format,
                                                  const QString& name)
{
    return {format, name,
            to_qstringlist(rdkit_extensions::get_rxn_extensions(format))};
}

} // unnamed namespace

// we define get_standard_formats and get_reaction_formats as functions rather
// than const values because they make use of const values defined in
// rdkit_extensions/file_format.cpp, and we can't guarantee that the file_format
// values will be defined before any values here get defined
FormatList<Format> get_standard_formats()
{
    return {
        // SKETCH-1453: Forbid MDL_MOLV2000 on export; potential stereo
        // ambiguities
        mol_data(Format::MDL_MOLV3000, "MDL SD V3000"),
        mol_data(Format::MAESTRO, "Maestro"),
        mol_data(Format::SMILES, "SMILES"),
        mol_data(Format::EXTENDED_SMILES, "Extended SMILES"),
        mol_data(Format::SMARTS, "SMARTS"),
        mol_data(Format::EXTENDED_SMARTS, "Extended SMARTS"),
        mol_data(Format::INCHI, "InChI"),
        mol_data(Format::INCHI_KEY, "InChIKey"),
        mol_data(Format::PDB, "PDB"),
        mol_data(Format::XYZ, "XYZ"),
    };
}

FormatList<Format> get_reaction_formats()
{
    return {
        rxn_data(Format::MDL_MOLV3000, "MDL RXN V3000"),
        rxn_data(Format::SMILES, "Reaction SMILES"),
        rxn_data(Format::SMARTS, "Reaction SMARTS"),
    };
}

// We define get_image_formats as a function for consistency with
// get_standard_formats and get_reaction_formats, even though it doesn't depend
// on any other values
FormatList<ImageFormat> get_image_formats()
{
    return {
        {ImageFormat::PNG, "PNG", {".png"}},
        {ImageFormat::SVG, "SVG", {".svg"}},
    };
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

QString get_import_name_filters()
{
    QStringList filters;
    for (const auto& format_list :
         {get_standard_formats(), get_reaction_formats()}) {
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
