/* -------------------------------------------------------------------------
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include "schrodinger/rdkit_extensions/file_format.h"

#include <fstream>
#include <unordered_map>

#include <boost/algorithm/string.hpp>
#include <fmt/format.h>

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
    {Format::SMILES, {".smi", ".smiles", ".smigz", ".smi.gz"}},
    {Format::EXTENDED_SMILES, {".cxsmi", ".cxsmiles"}},
    {Format::SMARTS, {}},
    {Format::EXTENDED_SMARTS, {}},
    {Format::MDL_MOLV2000, {}},
    {Format::MDL_MOLV3000,
     {".sdf", ".sd", ".mol", ".mdl", ".sdf.gz", ".sd.gz", ".mol.gz", ".sdfgz"}},
    {Format::MAESTRO, {".mae", ".maegz", ".mae.gz", ".maezst", ".mae.zst"}},
    {Format::INCHI, {".inchi"}},
    {Format::INCHI_KEY, {}},
    {Format::PDB, {".pdb", ".ent", ".pdb.gz", ".ent.gz", ".pdbgz", ".entgz"}},
    {Format::MOL2, {".mol2"}},
    {Format::XYZ, {".xyz"}},
    {Format::MRV, {".mrv"}},
    {Format::CDXML, {".cdxml"}},
};
// clang-format off
const std::unordered_map<Format, std::vector<std::string>> RXN_FORMAT_EXTS = {
    {Format::RDMOL_BINARY_BASE64, {}},
    {Format::SMILES, {".rsmi"}},
    {Format::EXTENDED_SMILES, {}},
    {Format::SMARTS, {}},
    {Format::EXTENDED_SMARTS, {}},
    {Format::MDL_MOLV2000, {}},
    {Format::MDL_MOLV3000, {".rxn"}},
};
// clang-format on
const std::unordered_map<Format, std::vector<std::string>> SEQ_FORMAT_EXTS = {
    {Format::RDMOL_BINARY_BASE64, {}},
    {Format::HELM, {".helm"}},
    {Format::FASTA_PEPTIDE, {}},
    {Format::FASTA_DNA, {}},
    {Format::FASTA_RNA, {}},
    {Format::FASTA, {".fasta", ".fas", ".fa", ".fst", ".seq"}},
};

std::vector<Format> get_roundtrip_keys(
    const std::unordered_map<Format, std::vector<std::string>>& map)
{
    std::vector<Format> keys;
    for (const auto& [key, _] : map) {
        // Exclude any formats that are not round-trippable
        if (key == Format::INCHI_KEY || key == Format::MOL2 ||
            key == Format::CDXML) {
            continue;
        }
        keys.push_back(key);
    }
    return keys;
}

} // unnamed namespace

// Vectors of round-trippable formats, primarily used in testing
const std::vector<Format> MOL_FORMATS = get_roundtrip_keys(MOL_FORMAT_EXTS);
const std::vector<Format> RXN_FORMATS = get_roundtrip_keys(RXN_FORMAT_EXTS);
const std::vector<Format> SEQ_FORMATS = get_roundtrip_keys(SEQ_FORMAT_EXTS);

std::vector<std::string> get_mol_extensions(const Format format)
{
    return MOL_FORMAT_EXTS.at(format);
}

std::vector<std::string> get_rxn_extensions(const Format format)
{
    return RXN_FORMAT_EXTS.at(format);
}

std::vector<std::string> get_seq_extensions(const Format format)
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

namespace
{
[[nodiscard]] CompressionType get_compression_type_impl(std::istream& is)
{
    constexpr unsigned char GZIP_MAGIC[2] = {0x1f, 0x8b};
    constexpr unsigned char ZSTD_MAGIC[4] = {0x28, 0xb5, 0x2f, 0xfd};

    // Read the first two bytes from the stream to determine the compression
    // type
    char magic[2];
    is.read(magic, 2);
    if (!is.good() || is.gcount() != 2) {
        return CompressionType::UNKNOWN;
    }

    auto u_magic = reinterpret_cast<unsigned char*>(magic);
    if (u_magic[0] == GZIP_MAGIC[0] && u_magic[1] == GZIP_MAGIC[1]) {
        return CompressionType::GZIP;
    } else if (u_magic[0] == ZSTD_MAGIC[0] && u_magic[1] == ZSTD_MAGIC[1]) {
        // zstd magic number is 4 bytes long, so read another 2 bytes
        is.read(magic, 2);
        if (is.gcount() == 2 && u_magic[0] == ZSTD_MAGIC[2] &&
            u_magic[1] == ZSTD_MAGIC[3]) {
            return CompressionType::ZSTD;
        }
    }
    return CompressionType::UNKNOWN;
}
} // namespace

CompressionType get_compression_type(const boost::filesystem::path& filename)
{
    std::ifstream file(filename.string(), std::ios::in | std::ios::binary);
    if (!file.good()) {
        throw fmt::system_error(errno, "Unable to open file {}",
                                filename.string());
    }

    return get_compression_type_impl(file);
}

CompressionType get_compression_type(std::istringstream& sstream)
{
    // save this so we can restore position on exit
    const auto prev_pos = sstream.tellg();

    // seek to the beginning of the stream
    sstream.seekg(0);
    if (!sstream) {
        throw fmt::system_error(errno, "Bad string stream in `{}`",
                                __FUNCTION__);
    }

    auto compression_type = get_compression_type_impl(sstream);
    // restore previous position
    sstream.seekg(prev_pos);

    return compression_type;
}

CompressionType
get_compression_type_from_ext(const boost::filesystem::path& filename)
{
    if (boost::algorithm::iends_with(filename.string(), "gz")) {
        return CompressionType::GZIP;
    } else if (boost::algorithm::iends_with(filename.string(), ".zst")) {
        return CompressionType::ZSTD;
    }
    return CompressionType::UNKNOWN;
}

std::ostream& operator<<(std::ostream& os,
                         const CompressionType& compression_type)
{
    switch (compression_type) {
        case CompressionType::GZIP:
            return os << "CompressionType::GZIP";
        case CompressionType::ZSTD:
            return os << "CompressionType::ZSTD";
        case CompressionType::UNKNOWN:
            return os << "CompressionType::UNKNOWN";
    }
    return os;
}

} // namespace rdkit_extensions
} // namespace schrodinger
