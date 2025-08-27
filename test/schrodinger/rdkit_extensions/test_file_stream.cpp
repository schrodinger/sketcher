/* -------------------------------------------------------------------------
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_MODULE test_file_stream

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <cstdio>
#include <sstream>
#include <ios>

#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/rdkit_extensions/file_stream.h"
#include "schrodinger/test/checkexceptionmsg.h"
#include "test_common.h"

using namespace schrodinger::rdkit_extensions;
namespace bdata = boost::unit_test::data;

BOOST_AUTO_TEST_CASE(TestReadingUncompressedFiles)
{
    auto fname = testfile_path("methane.sdf");
    maybe_compressed_istream fstream(fname);

    BOOST_TEST(fstream.good()); // should be readable
    BOOST_TEST(!fstream.is_compressed());

    constexpr size_t read_size = 100;
    std::string buffer(read_size, '\0');
    fstream.read(&buffer[0], read_size);

    BOOST_TEST(fstream.gcount() == read_size);
    BOOST_TEST(fstream.tellg() == read_size);
    BOOST_TEST(buffer != std::string(read_size, '\0'));

    BOOST_TEST(fstream.good()); // should still be readable
}

BOOST_DATA_TEST_CASE(TestReadingCompressedFiles,
                     boost::unit_test::data::make(std::vector<std::string>{
                         "methane.smigz", "methane.sdfgz", "methane.maegz",
                         "methane.mae.zst"}),
                     testfile)
{
    auto fname = testfile_path(testfile);
    maybe_compressed_istream fstream(fname);

    BOOST_TEST(fstream.good()); // should be readable
    BOOST_TEST(fstream.is_compressed());

    constexpr size_t read_size = 2;
    std::string buffer(read_size, '\0');
    fstream.read(&buffer[0], read_size);

    // can't test tellg() since output passes through multiple buffer layers
    BOOST_TEST(fstream.gcount() == read_size);
    BOOST_TEST(buffer != std::string(read_size, '\0'));

    BOOST_TEST(fstream.good()); // should still be readable
}

BOOST_DATA_TEST_CASE(
    TestWritingInterface,
    boost::unit_test::data::make(std::vector<std::string>{
        ".mae",
        ".mae.gz",
        ".sdf",
        ".sdfgz",
        ".mae.zst",
    }) * boost::unit_test::data::make(std::vector<std::ios::openmode>{
             std::ios::app, std::ios::trunc}),
    test_ext, write_mode)
{
    auto tmpfile = boost::filesystem::temp_directory_path() /
                   boost::filesystem::unique_path();
    auto testfile = tmpfile.string() + test_ext;

    maybe_compressed_ostream fstream(testfile, write_mode);
    BOOST_TEST(fstream.good()); // should be writable

    fstream << "hello there\n";
    BOOST_TEST(fstream.good()); // should still be writable
}

BOOST_DATA_TEST_CASE(TestReadingFromStringInput,
                     bdata::make(std::vector<std::string>{
                         "methane.sdf", "methane.smigz", "methane.sdfgz",
                         "methane.maegz", "methane.mae.zst"}),
                     testfile)
{
    auto fname = testfile_path(testfile);
    std::ifstream is(fname, std::ios::binary);
    std::string text(std::istreambuf_iterator<char>(is), {});

    // open file stream with text
    auto compression_type = get_compression_type(fname);
    maybe_compressed_istream fstream(text, compression_type);
    BOOST_TEST(fstream.good()); // should be readable
    auto is_compressed = compression_type != CompressionType::UNKNOWN;
    BOOST_TEST(fstream.is_compressed() == is_compressed);

    constexpr size_t read_size = 2;
    std::string buffer(read_size, '\0');
    fstream.read(&buffer[0], read_size);

    // can't test tellg() since output passes through multiple buffer layers
    BOOST_TEST(fstream.gcount() == read_size);
    BOOST_TEST(buffer != std::string(read_size, '\0'));

    BOOST_TEST(fstream.good()); // should still be readable
}

BOOST_AUTO_TEST_CASE(TestInputStreamFailsForNonexistentFile)
{
    std::string fname = std::string(__FILE__) + ".bad_ext";
    TEST_CHECK_EXCEPTION_MSG_SUBSTR((maybe_compressed_istream(fname)),
                                    std::system_error, "Error opening");
}

BOOST_AUTO_TEST_CASE(TestOutputStreamFailsForNonexistentFolder)
{
    std::string fname = std::string(__FILE__) + "/dummy.txt";
    TEST_CHECK_EXCEPTION_MSG_SUBSTR((maybe_compressed_ostream(fname)),
                                    std::system_error, "Error opening");
}

BOOST_AUTO_TEST_CASE(TestOutputStreamFailsForBadOstream)
{
    std::ostringstream os;
    os.setstate(std::ios_base::failbit);

    TEST_CHECK_EXCEPTION_MSG_SUBSTR(
        maybe_compressed_ostream(os, CompressionType::UNKNOWN),
        std::runtime_error, "Bad output stream");
}

BOOST_AUTO_TEST_CASE(TestGetCompressedString)
{
    auto fname = testfile_path("methane.mae");
    std::ifstream is(fname);
    std::string buffer(std::istreambuf_iterator<char>(is), {});

    // NOTE: Compression type should be correct
    auto gz_string = get_compressed_string(buffer, CompressionType::GZIP);
    std::istringstream ss1(gz_string);
    BOOST_TEST(get_compression_type(ss1) == CompressionType::GZIP);

    auto zstd_string = get_compressed_string(buffer, CompressionType::ZSTD);
    std::istringstream ss2(zstd_string);
    BOOST_TEST(get_compression_type(ss2) == CompressionType::ZSTD);
}

BOOST_DATA_TEST_CASE(TestGetDecompressedString,
                     bdata::make(std::vector<std::string>{
                         "methane.mae", "methane.mae.zst", "methane.maegz"}) ^
                         bdata::make(std::vector<CompressionType>{
                             CompressionType::UNKNOWN, CompressionType::ZSTD,
                             CompressionType::GZIP}),
                     testfile, compression_type)
{
    auto fname = testfile_path(testfile);
    std::ifstream is(fname, std::ios::binary);
    std::string buffer(std::istreambuf_iterator<char>(is), {});

    // MAESTRO formatted file should have a title
    auto decompressed_string =
        get_decompressed_string(buffer, compression_type);
    BOOST_TEST(decompressed_string.find("s_m_title") != std::string::npos);

    decompressed_string = get_decompressed_string(buffer);
    BOOST_TEST(decompressed_string.find("s_m_title") != std::string::npos);
}
