#pragma once

/**
 * Open a named file for reading or writing. If the filename ends in .gz,
 * then it will be opened at gzipped compressed
 *
 * @copyright Copyright Schrödinger (c) 2023
 *
 */

#include <fstream>
#include <optional>
#include <sstream>

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"
#include "schrodinger/rdkit_extensions/file_format.h"

namespace schrodinger
{
namespace rdkit_extensions
{

// Open a filename for writing. "gz" and ".zst" suffixed files will be
// compressed
//
// @throws std::system_error if directory does not exist or if output stream
// is in a "bad" or "fail" state
class RDKIT_EXTENSIONS_API maybe_compressed_ostream : public std::ostream,
                                                      public boost::noncopyable
{
  public:
    explicit maybe_compressed_ostream(
        const boost::filesystem::path& filename,
        std::ios::openmode mode = std::ios::trunc);

    explicit maybe_compressed_ostream(std::ostream& output_stream,
                                      CompressionType compression_type);

  private:
    // this is needed if we're writing to disk
    std::optional<std::ofstream> m_sdgr_ofstream;
    boost::iostreams::filtering_ostreambuf m_sdgr_buffer;

    void initialize_ostream(std::ostream& os,
                            const CompressionType& compression_type);
};

// an interface to opening files for reading. This allows file reading logic for
// compressed and uncompressed inputs and will be a good place to implement the
// read progress logic for structure readers
//
// @throws std::system_error if file does not exist or if we can't create a
// string stream from the input buffer.
class RDKIT_EXTENSIONS_API maybe_compressed_istream : public std::istream,
                                                      public boost::noncopyable
{
  public:
    explicit maybe_compressed_istream(const boost::filesystem::path& filename);
    explicit maybe_compressed_istream(const std::string& data,
                                      CompressionType compression_type);

    bool is_compressed() const;

    //
    // returns the current position in the input stream.
    //
    // NOTE: compressed inputs are decompressed in chunks, so this allows us to
    // get the current position in the compressed input. The returned
    // position could be the same after reading from the input stream multiple
    // times.
    std::istream::pos_type tellg();

  private:
    // we always need a data stream
    std::optional<std::ifstream> m_sdgr_ifstream;
    std::optional<std::istringstream> m_sdgr_istringstream;

    // only compressed inputs need a filtering stream buffer
    std::optional<boost::iostreams::filtering_istreambuf> m_sdgr_buffer;
    bool m_sdgr_is_compressed = false;

    void initialize_istream(std::istream& is,
                            const CompressionType& compression_type);
};

[[nodiscard]] RDKIT_EXTENSIONS_API std::string
get_compressed_string(const std::string& data,
                      CompressionType compression_type);

[[nodiscard]] RDKIT_EXTENSIONS_API std::string
get_decompressed_string(const std::string& data);

[[nodiscard]] RDKIT_EXTENSIONS_API std::string
get_decompressed_string(const std::string& data,
                        CompressionType compression_type);
} // namespace rdkit_extensions
} // namespace schrodinger
