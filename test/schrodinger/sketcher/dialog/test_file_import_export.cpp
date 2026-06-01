
#define BOOST_TEST_MODULE test_file_import_export

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/test/checkexceptionmsg.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);

BOOST_AUTO_TEST_CASE(test_get_import_export_formats)
{
    BOOST_TEST(get_import_formats().size() == 12);
    // 11 atomistic + HELM + FASTA = 13
    BOOST_TEST(get_standard_export_formats().size() == 13);
    BOOST_TEST(get_reaction_export_formats().size() == 5);
    BOOST_TEST(get_image_export_formats().size() == 2);
}

BOOST_AUTO_TEST_CASE(test_get_standard_export_formats_includes_sequence_formats)
{
    auto formats = get_standard_export_formats();
    auto has_format = [&formats](Format target) {
        return std::ranges::any_of(formats, [target](const auto& e) {
            return std::get<0>(e) == target;
        });
    };
    BOOST_TEST(has_format(Format::HELM));
    BOOST_TEST(has_format(Format::FASTA));
}
