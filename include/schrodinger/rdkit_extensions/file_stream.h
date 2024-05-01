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

#include <boost/filesystem.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include "schrodinger/rdkit_extensions/definitions.h"

namespace schrodinger
{
namespace rdkit_extensions
{

// Open a filename for writing. "gz" and ".zst" suffixed files will be
// compressed
class RDKIT_EXTENSIONS_API maybe_compressed_ostream : public std::ostream,
                                                      public boost::noncopyable
{
  public:
    explicit maybe_compressed_ostream(
        const boost::filesystem::path& filename,
        std::ios::openmode mode = std::ios::trunc);

  private:
    boost::iostreams::filtering_ostreambuf m_sdgr_buffer;
};

// an interface to opening files for reading. This allows file reading logic for
// compressed and uncompressed inputs and will be a good place to implement the
// read progress logic for structure readers
class RDKIT_EXTENSIONS_API maybe_compressed_istream : public std::istream,
                                                      public boost::noncopyable
{
  public:
    explicit maybe_compressed_istream(const boost::filesystem::path& filename);

    bool is_compressed() const;

    //
    // returns the current file position in the input file stream.
    //
    // NOTE: compressed inputs are decompressed in chunks, so this allows us to
    // get the current file position in the compressed input. The returned
    // position could be the same after reading from the file multiple times.
    std::istream::pos_type tellg();

  private:
    // we always need a file stream
    std::ifstream m_sdgr_ifstream;
    // only compressed inputs need a filtering stream buffer
    std::optional<boost::iostreams::filtering_istreambuf> m_sdgr_buffer;
    bool m_sdgr_is_compressed = false;
};
} // namespace rdkit_extensions
} // namespace schrodinger
