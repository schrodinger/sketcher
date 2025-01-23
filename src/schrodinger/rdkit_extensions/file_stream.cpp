#include "schrodinger/rdkit_extensions/file_stream.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zstd.hpp>
#include <cstring>
#include <fmt/format.h>

#include "schrodinger/rdkit_extensions/file_format.h"

namespace schrodinger
{
namespace rdkit_extensions
{

maybe_compressed_ostream::maybe_compressed_ostream(
    const boost::filesystem::path& filename, std::ios::openmode mode) :
    std::ostream(nullptr)
{
    auto fpath = filename.string();
    m_sdgr_ofstream.emplace(fpath,
                            mode | std::ios_base::out | std::ios::binary);
    if (!*m_sdgr_ofstream) {
        throw fmt::system_error(errno, "Error opening {}", fpath);
    }

    initialize_ostream(*m_sdgr_ofstream,
                       get_compression_type_from_ext(filename));
}

maybe_compressed_ostream::maybe_compressed_ostream(
    std::ostream& output_stream, CompressionType compression_type) :
    std::ostream(nullptr)
{
    if (!output_stream) {
        throw std::runtime_error(
            "Bad output stream in `maybe_compressed_ostream`");
    }

    initialize_ostream(output_stream, compression_type);
}

void maybe_compressed_ostream::initialize_ostream(
    std::ostream& output_stream, const CompressionType& compression_type)
{
    if (compression_type == CompressionType::GZIP) {
        m_sdgr_buffer.push(boost::iostreams::gzip_compressor());
    } else if (compression_type == CompressionType::ZSTD) {
        m_sdgr_buffer.push(boost::iostreams::zstd_compressor());
    }

    m_sdgr_buffer.push(boost::ref(output_stream));
    this->init(&m_sdgr_buffer);

    // NOTE: When we encounter a bad write stream state, we want to throw an
    // error instead of just setting the fail or bad bit. This is important for
    // cases like the one reported in SHARED-11255 since the write stream could
    // write a corrupted output file and fail silently without any message to
    // the user.
    constexpr auto error_states = std::ios::failbit | std::ios::badbit;
    // make the underlying stream throw if the filtering stream fails to catch
    // the relevant error.
    output_stream.exceptions(error_states);
    // the filtering stream should throw on error
    this->exceptions(error_states);
}

maybe_compressed_istream::maybe_compressed_istream(
    const boost::filesystem::path& filename) :
    std::istream(nullptr)
{
    // NOTE: We should always open in binary to make file seeking consistent
    // across different platforms. See SHARED-10715
    auto fpath = filename.string();
    m_sdgr_ifstream.emplace(fpath, std::ios::in | std::ios::binary);
    if (!*m_sdgr_ifstream) {
        throw fmt::system_error(errno, "Error opening {}", fpath);
    }

    auto compression_type = get_compression_type(filename);
    initialize_istream(*m_sdgr_ifstream, compression_type);
}

maybe_compressed_istream::maybe_compressed_istream(
    const std::string& data, CompressionType compression_type) :
    std::istream(nullptr)
{
    m_sdgr_istringstream.emplace(data);
    if (!*m_sdgr_istringstream) {
        throw fmt::system_error(errno, "Error creating string stream");
    }

    initialize_istream(*m_sdgr_istringstream, compression_type);
}

bool maybe_compressed_istream::is_compressed() const
{
    return m_sdgr_is_compressed;
}

std::istream::pos_type maybe_compressed_istream::tellg()
{
    // assume one of these is not empty
    return m_sdgr_ifstream ? m_sdgr_ifstream->tellg()
                           : m_sdgr_istringstream->tellg();
}

void maybe_compressed_istream::initialize_istream(
    std::istream& is, const CompressionType& compression_type)
{
    m_sdgr_is_compressed = compression_type != CompressionType::UNKNOWN;
    // if we don't do this, i.e. if we use the filtering buffer for reading
    // uncompressed inputs, we won't have a useful tellg function because
    // we'll have multiple buffers between the file stream and the output
    if (!m_sdgr_is_compressed) {
        this->init(is.rdbuf());
        return;
    }

    // From Dan's benchmarks of the MaeReader
    constexpr auto BUFFER_SIZE_FOR_COMPRESSED = 4098 * 8;

    m_sdgr_buffer.emplace();
    if (compression_type == CompressionType::GZIP) {
        m_sdgr_buffer->push(boost::iostreams::gzip_decompressor(),
                            BUFFER_SIZE_FOR_COMPRESSED);
    } else if (compression_type == CompressionType::ZSTD) {
        m_sdgr_buffer->push(boost::iostreams::zstd_decompressor(),
                            BUFFER_SIZE_FOR_COMPRESSED);
    } else {
        // this is currently never hit, but we should be aware if we extend
        // the supported compression types
        throw std::invalid_argument("Unsupported compression type");
    }

    m_sdgr_buffer->push(boost::ref(is), BUFFER_SIZE_FOR_COMPRESSED);
    this->init(&*m_sdgr_buffer);
}

std::string get_compressed_string(const std::string& data,
                                  CompressionType compression_type)
{
    if (compression_type == CompressionType::UNKNOWN) {
        return data;
    } else {
        std::stringstream ss;
        // use maybe_compressed_ostream as RAII
        {
            maybe_compressed_ostream os(ss, compression_type);
            os.write(data.data(), data.size());
        }
        return ss.str();
    }
}

std::string get_decompressed_string(const std::string& data,
                                    CompressionType compression_type)
{
    if (compression_type == CompressionType::UNKNOWN) {
        return data;
    } else {
        maybe_compressed_istream is(data, compression_type);
        return std::string(std::istreambuf_iterator<char>(is), {});
    }
}

std::string get_decompressed_string(const std::string& data)
{
    std::istringstream sstream(data);
    auto compression_type = get_compression_type(sstream);
    return get_decompressed_string(data, compression_type);
}

} // namespace rdkit_extensions
} // namespace schrodinger
