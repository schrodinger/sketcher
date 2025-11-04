/* -------------------------------------------------------------------------
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_MODULE test_file_format

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>
#include <fmt/format.h>
#include <fstream>

#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/test/checkexceptionmsg.h"
#include "test_common.h"

using namespace schrodinger::rdkit_extensions;

BOOST_TEST_DONT_PRINT_LOG_VALUE(CompressionType);
BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);

/**
 * These test roundtrip all extensions through the public APIs to make sure
 * that all extensions are uniquely associated with a single format.
 */
BOOST_DATA_TEST_CASE(test_get_mol_extensions,
                     boost::unit_test::data::make(MOL_FORMATS))
{
    for (const auto& file_extension : get_mol_extensions(sample)) {
        BOOST_TEST(get_file_format(file_extension) == sample);
    }
}

BOOST_DATA_TEST_CASE(test_get_rxn_extensions,
                     boost::unit_test::data::make(RXN_FORMATS))
{
    for (const auto& file_extension : get_rxn_extensions(sample)) {
        BOOST_TEST(get_file_format(file_extension) == sample);
    }
}

BOOST_DATA_TEST_CASE(test_get_seq_extensions,
                     boost::unit_test::data::make(SEQ_FORMATS))
{
    for (const auto& file_extension : get_seq_extensions(sample)) {
        BOOST_TEST(get_file_format(file_extension) == sample);
    }
}

BOOST_AUTO_TEST_CASE(test_unsupported_extensions)
{
    // Incomplete and unknown extensions throw
    for (const auto& filename :
         {"", ".", ".foo", "pdb", "test.cx", "test.inchikey"}) {
        TEST_CHECK_EXCEPTION_MSG_SUBSTR(get_file_format(filename),
                                        std::invalid_argument,
                                        "Unsupported file extension");
    }
}

BOOST_DATA_TEST_CASE(
    test_compressed_formats,
    boost::unit_test::data::make({Format::SMILES, Format::MDL_MOLV3000,
                                  Format::MAESTRO, Format::PDB}),
    format)
{
    const auto extension = get_mol_extensions(format).front();

    auto filename = fmt::format("test{}", extension);
    BOOST_TEST(get_file_format(filename) == format);

    filename = fmt::format("test{}gz", extension);
    BOOST_TEST(get_file_format(filename) == format);

    filename = fmt::format("test{}.gz", extension);
    BOOST_TEST(get_file_format(filename) == format);
}

BOOST_DATA_TEST_CASE(
    TestGetCompressionType,
    boost::unit_test::data::make(std::vector<std::string>{
        ".mae", ".mol2", ".pdb", ".sdf", ".smi", ".smigz", ".sdfgz", ".maegz",
        ".mae.zst"}) ^
        boost::unit_test::data::make(std::vector<CompressionType>{
            CompressionType::UNKNOWN,
            CompressionType::UNKNOWN,
            CompressionType::UNKNOWN,
            CompressionType::UNKNOWN,
            CompressionType::UNKNOWN,
            CompressionType::GZIP,
            CompressionType::GZIP,
            CompressionType::GZIP,
            CompressionType::ZSTD,
        }),
    extension, expected_compression_type)
{

    auto fname = testfile_path("methane" + extension);
    BOOST_TEST(get_compression_type(fname) == expected_compression_type);
    BOOST_TEST(get_compression_type_from_ext(fname) ==
               expected_compression_type);

    // we should be able to determine the compression type from a
    // string stream
    std::ifstream is(fname);
    std::string buffer(std::istreambuf_iterator<char>(is), {});
    std::istringstream ss(buffer);
    BOOST_TEST(get_compression_type(ss) == expected_compression_type);
}
