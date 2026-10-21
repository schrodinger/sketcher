
#define BOOST_TEST_MODULE test_file_import_export

#include <boost/test/unit_test.hpp>

#include <QStringList>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/test/checkexceptionmsg.h"

using namespace schrodinger::sketcher;
using schrodinger::rdkit_extensions::Format;

BOOST_TEST_DONT_PRINT_LOG_VALUE(Format);

BOOST_AUTO_TEST_CASE(test_get_import_export_formats)
{
    // 10 atomistic (SMARTS variants have no extensions) + 2 reaction + HELM
    BOOST_TEST(get_import_formats().size() == 13);
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

BOOST_AUTO_TEST_CASE(test_get_import_formats_includes_helm)
{
    // SKETCH-2823: HELM files should be selectable in the import dialog
    auto formats = get_import_formats();
    BOOST_TEST(std::ranges::any_of(
        formats, [](const auto& e) { return std::get<0>(e) == Format::HELM; }));
}

BOOST_AUTO_TEST_CASE(test_get_import_name_filter)
{
    auto filters = get_import_name_filter().split(";;");
#ifdef __EMSCRIPTEN__
    // SKETCH-2823: Qt drops wildcard filters on WASM, so any other filter
    // would prevent the browser from showing all files
    BOOST_TEST((filters == QStringList{"All files (*)"}));
#else
    // SKETCH-2823: All supported formats, then each format
    BOOST_TEST(filters.size() == get_import_formats().size() + 1);
    auto all_supported = filters.front();
    BOOST_TEST(all_supported.startsWith("All supported formats ("));
    BOOST_TEST(!all_supported.contains("*.*"));
    for (const auto& [_, label, extensions] : get_import_formats()) {
        for (const auto& ext : extensions) {
            BOOST_TEST(all_supported.contains(QString::fromStdString(ext)));
        }
    }
#endif
}
