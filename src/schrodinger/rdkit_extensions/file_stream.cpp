#include "schrodinger/rdkit_extensions/file_stream.h"

#include <boost/algorithm/string/predicate.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filter/zstd.hpp>

#include "schrodinger/rdkit_extensions/file_format.h"

namespace schrodinger
{
namespace rdkit_extensions
{

maybe_compressed_ostream::maybe_compressed_ostream(
    const boost::filesystem::path& filename, std::ios::openmode mode) :
    std::ostream(nullptr),
    m_sdgr_buffer()
{
    if (boost::algorithm::iends_with(filename.string(), "gz")) {
        m_sdgr_buffer.push(boost::iostreams::gzip_compressor());
    } else if (boost::algorithm::iends_with(filename.string(), ".zst")) {
        m_sdgr_buffer.push(boost::iostreams::zstd_compressor());
    }

    // turn on the required flags
    mode |= (std::ios_base::out | std::ios::binary);
    m_sdgr_buffer.push(boost::iostreams::file_descriptor_sink(filename, mode));
    this->init(&m_sdgr_buffer);
}

maybe_compressed_istream::maybe_compressed_istream(
    const boost::filesystem::path& filename) :
    std::istream(nullptr),
    // NOTE: We should always open in binary to make file seeking consistent
    // across different platforms. See SHARED-10715
    m_sdgr_ifstream(filename.string(), std::ios::in | std::ios::binary)
{
    const auto compression_type = get_compression_type(filename);
    // if we don't do this, i.e. if we use the filtering buffer for reading
    // uncompressed inputs, we won't have a useful tellg function because
    // we'll have multiple buffers between the file stream and the output
    if (compression_type == CompressionType::UNKNOWN) {
        this->init(m_sdgr_ifstream.rdbuf());
        m_sdgr_is_compressed = false;
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

    m_sdgr_buffer->push(boost::ref(m_sdgr_ifstream),
                        BUFFER_SIZE_FOR_COMPRESSED);
    this->init(&*m_sdgr_buffer);
    m_sdgr_is_compressed = true;
}

bool maybe_compressed_istream::is_compressed() const
{
    return m_sdgr_is_compressed;
}

std::istream::pos_type maybe_compressed_istream::tellg()
{
    return m_sdgr_ifstream.tellg();
}

} // namespace rdkit_extensions
} // namespace schrodinger
