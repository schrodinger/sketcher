/* -------------------------------------------------------------------------
 * Tests schrodinger::sketcher:: rendering APIs
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sketcher_image_generation

#include <boost/filesystem.hpp>
#include <boost/test/framework.hpp>
#include <boost/test/unit_test.hpp>

#include <GraphMol/GraphMol.h>

#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "test_common.h"

using boost::unit_test::framework::current_test_case;
using namespace schrodinger::sketcher;
using namespace schrodinger;

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

BOOST_AUTO_TEST_CASE(test_image_gen_APIs)
{
    auto rdmol = rdkit_extensions::text_to_rdmol("C1=CC=CC=C1");

    RenderOptions opts;
    opts.width_height = {400, 400};
    opts.background_color = QColor("#d4e6f1");

    auto qpict = get_qpicture(*rdmol, opts);
    BOOST_TEST(qpict.size() > 0);

    auto image = get_qimage(*rdmol, opts);
    BOOST_TEST(image.sizeInBytes() > 0);

    auto bytes = get_image_bytes(*rdmol, ImageFormat::PNG, opts);
    BOOST_TEST(bytes.size() > 0);

    const auto tmp_file = "tmp_for_" + current_test_case().p_name.get();
    for (const auto& ext : {".png", ".svg"}) {
        auto outfile = tmp_file + ext;
        save_image_file(*rdmol, outfile, opts);
        BOOST_TEST(boost::filesystem::exists(outfile));
    }
    BOOST_REQUIRE_THROW(save_image_file(*rdmol, tmp_file + ".invalid", opts),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_highlighting)
{
    auto rdmol = rdkit_extensions::text_to_rdmol("C1=CC=CC=C1");

    RenderOptions opts;
    opts.width_height = {400, 400};
    opts.background_color = QColor(255, 255, 255);
    QHash<int, QColor> highlights;
    QColor red(255, 0, 0);
    highlights[0] = red;
    highlights[1] = red;
    highlights[2] = red;
    QColor purple(255, 0, 25);
    highlights[3] = purple;
    highlights[4] = purple;
    opts.rdatom_index_to_highlight_color = highlights;

    QHash<int, QColor> highlights_2;
    QColor yellow(255, 255, 0);
    highlights[0] = yellow;
    highlights[1] = yellow;

    opts.rdbond_index_to_highlight_color = highlights_2;

    auto qimage = get_qimage(*rdmol, opts);
    BOOST_TEST(!qimage.allGray());
}
