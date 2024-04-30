/* -------------------------------------------------------------------------
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include "schrodinger/rdkit_extensions/file_format.h"

#include <iostream>
#include <unordered_map>

#include <boost/algorithm/string.hpp>

namespace schrodinger
{
namespace rdkit_extensions
{

namespace
{
/**
 * @internal
 * These maps are set up such that keys indicate all formats that specific to
 * what is supported when reading/writing either as RDKit mols or chemical
 * reactions. If a Format is omitted from the list, that means there is
 * currently no way to represent either a mol or reaction in that format. If
 * the format indexes an empty list, that means there is either no typical file
 * extension associated with that format, or there is a preferred format for
 * writing that format type. The latter is applicable to MDL_MOLV2000 (where
 * MDL_MOLV3000 is preferred) and the various specific FASTA sub-formats.
 *
 * https://openbabel.org/docs/current/FileFormats/Common_cheminformatics_Formats.html
 * https://search.r-project.org/CRAN/refmans/seqinr/html/read.fasta.html
 */
const std::unordered_map<Format, std::vector<std::string>> MOL_FORMAT_EXTS = {
    {Format::RDMOL_BINARY_BASE64, {}},
    {Format::SMILES, {".smi", ".smiles"}},
    {Format::EXTENDED_SMILES, {".cxsmi", ".cxsmiles"}},
    {Format::SMARTS, {}},
    {Format::MDL_MOLV2000, {}},
    {Format::MDL_MOLV3000, {".sdf", ".sd", ".mol", ".mdl"}},
    {Format::MAESTRO, {".mae"}},
    {Format::INCHI, {".inchi"}},
    {Format::INCHI_KEY, {}},
    {Format::PDB, {".pdb", ".ent"}},
    {Format::XYZ, {".xyz"}},
};

const std::unordered_map<Format, std::vector<std::string>> RXN_FORMAT_EXTS = {
    {Format::RDMOL_BINARY_BASE64, {}},
    {Format::SMILES, {".rsmi"}},
    {Format::SMARTS, {}},
    {Format::MDL_MOLV2000, {}},
    {Format::MDL_MOLV3000, {".rxn"}},
};

const std::unordered_map<Format, std::vector<std::string>> SEQ_FORMAT_EXTS = {
    {Format::RDMOL_BINARY_BASE64, {}},
    {Format::HELM, {".helm"}},
    {Format::FASTA_PEPTIDE, {}},
    {Format::FASTA_DNA, {}},
    {Format::FASTA_RNA, {}},
    {Format::FASTA, {".fasta", ".fas", ".fa", ".fst", ".seq"}},
};

std::vector<Format>
get_keys(const std::unordered_map<Format, std::vector<std::string>>& map)
{
    std::vector<Format> keys;
    for (const auto& [key, _] : map) {
        keys.push_back(key);
    }
    return keys;
}

} // unnamed namespace

const std::vector<Format> MOL_FORMATS = get_keys(MOL_FORMAT_EXTS);
const std::vector<Format> RXN_FORMATS = get_keys(RXN_FORMAT_EXTS);
const std::vector<Format> SEQ_FORMATS = get_keys(SEQ_FORMAT_EXTS);

std::vector<std::string> get_mol_extensions(Format format)
{
    return MOL_FORMAT_EXTS.at(format);
}

std::vector<std::string> get_rxn_extensions(Format format)
{
    return RXN_FORMAT_EXTS.at(format);
}

std::vector<std::string> get_seq_extensions(Format format)
{
    return SEQ_FORMAT_EXTS.at(format);
}

Format get_file_format(const boost::filesystem::path& filename)
{
    auto extension = filename.extension().string();
    boost::to_lower(extension);

    if (extension == ".gz" || extension == ".zst") {
        extension = filename.stem().extension().string() + extension;
        boost::to_lower(extension);
    }

    for (const auto& format_to_extensions_map :
         {MOL_FORMAT_EXTS, RXN_FORMAT_EXTS, SEQ_FORMAT_EXTS}) {
        for (const auto& [format, exts] : format_to_extensions_map) {
            if (std::find(exts.begin(), exts.end(), extension) != exts.end()) {
                return format;
            }
        }
    }
    throw std::invalid_argument("Unsupported file extension: " +
                                filename.string());
}

} // namespace rdkit_extensions
} // namespace schrodinger
