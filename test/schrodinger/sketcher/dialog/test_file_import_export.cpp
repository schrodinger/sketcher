#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_file_import_export

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/test/checkexceptionmsg.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);

BOOST_AUTO_TEST_CASE(test_get_file_format)
{
    BOOST_TEST(get_file_format("test.smi") == Format::SMILES);
    BOOST_TEST(get_file_format("test.smiles") == Format::SMILES);
    BOOST_TEST(get_file_format("test.cxsmi") == Format::EXTENDED_SMILES);
    BOOST_TEST(get_file_format("test.cxsmiles") == Format::EXTENDED_SMILES);
    BOOST_TEST(get_file_format("test.sdf") == Format::MDL_MOLV3000);
    BOOST_TEST(get_file_format("test.sd") == Format::MDL_MOLV3000);
    BOOST_TEST(get_file_format("test.mol") == Format::MDL_MOLV3000);
    BOOST_TEST(get_file_format("test.mdl") == Format::MDL_MOLV3000);
    BOOST_TEST(get_file_format("test.mae") == Format::MAESTRO);
    BOOST_TEST(get_file_format("test.maegz") == Format::MAESTRO);
    BOOST_TEST(get_file_format("test.mae.gz") == Format::MAESTRO);
    BOOST_TEST(get_file_format("test.inchi") == Format::INCHI);
    BOOST_TEST(get_file_format("test.pdb") == Format::PDB);
    BOOST_TEST(get_file_format("test.ent") == Format::PDB);
    BOOST_TEST(get_file_format("test.rsmi") == Format::SMILES);
    BOOST_TEST(get_file_format("test.rxn") == Format::MDL_MOLV3000);

    // Incomplete and unknown extensions throw
    for (const auto& filename :
         {"test.cx", "test.sdf.gz", "test.inchikey", "test.with.dots.rxn"}) {
        TEST_CHECK_EXCEPTION_MSG_SUBSTR(get_file_format(filename),
                                        std::runtime_error,
                                        "Unknown file extension");
    }

    auto filters = get_import_name_filters();
    BOOST_TEST(filters.contains(";;"));
    BOOST_TEST(filters.split(";;").size() == 8);
    BOOST_TEST(filters.contains("*.sdf"));
    BOOST_TEST(filters.contains("*.mdl"));
}
