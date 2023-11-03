/* -------------------------------------------------------------------------
 * Tests schrodinger::sketcher:: rendering APIs
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE sketcher_image_generation

#include <GraphMol/GraphMol.h>
#include <QList>
#include <boost/filesystem.hpp>
#include <boost/test/framework.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/image_generation.h"
#include "test_common.h"

using boost::unit_test::framework::current_test_case;
using namespace schrodinger::sketcher;
using namespace schrodinger;

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

template <typename T> void image_gen_APIs(const T& input)
{
    RenderOptions opts;
    opts.width_height = {400, 400};
    opts.background_color = QColor("#d4e6f1");

    auto qpict = get_qpicture(input, opts);
    BOOST_TEST(qpict.size() > 0);

    auto image = get_qimage(input, opts);
    BOOST_TEST(image.sizeInBytes() > 0);

    auto bytes = get_image_bytes(input, ImageFormat::PNG, opts);
    BOOST_TEST(bytes.size() > 0);

    const auto tmp_file = "tmp_for_" + current_test_case().p_name.get();
    for (const auto& ext : {".png", ".svg"}) {
        auto outfile = tmp_file + ext;
        save_image_file(input, outfile, opts);
        BOOST_TEST(boost::filesystem::exists(outfile));
    }
    BOOST_REQUIRE_THROW(save_image_file(input, tmp_file + ".invalid", opts),
                        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(test_image_gen_APIs)
{
    auto mol_smiles = "C1=CC=CC=C1";
    auto rwmol = rdkit_extensions::to_rdkit(mol_smiles);
    image_gen_APIs(*rwmol);
    image_gen_APIs(RDKit::ROMol(*rwmol));
    image_gen_APIs(mol_smiles);
    image_gen_APIs(std::string(mol_smiles));

    auto rxn_smiles = "CC(=O)O.OCC>>CC(=O)OCC";
    image_gen_APIs(*rdkit_extensions::to_rdkit_reaction(rxn_smiles));
    image_gen_APIs(rxn_smiles);
    image_gen_APIs(std::string(rxn_smiles));
}

BOOST_AUTO_TEST_CASE(test_highlighting)
{
    auto rdmol = rdkit_extensions::to_rdkit("C1=CC=CC=C1");

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

BOOST_AUTO_TEST_CASE(test_get_best_image_scale)
{
    QList<RDKit::ROMol*> rdmols;
    BOOST_TEST(get_best_image_scale(rdmols) == -1);
    auto rdmol_small = rdkit_extensions::to_rdkit("c");
    auto rdmol_med = rdkit_extensions::to_rdkit("c1ccccc1");
    auto rdmol_large = rdkit_extensions::to_rdkit("c1nccc2n1ccc2");
    rdmols.append({rdmol_small.get(), rdmol_med.get(), rdmol_large.get()});
    qreal scale_all_three = get_best_image_scale(rdmols);
    BOOST_TEST(scale_all_three > 0);
    qreal scale_small = get_best_image_scale({rdmol_small.get()});
    BOOST_TEST(scale_small > 0);
    // removing larger molecules should increase the scale
    BOOST_TEST(scale_small > scale_all_three);
    qreal scale_large = get_best_image_scale({rdmol_large.get()});
    // using only the largest molecule should give the same scale
    BOOST_CHECK_CLOSE(scale_large, scale_all_three, 0.01);
    RenderOptions opts;
    opts.width_height = {100, 100};
    qreal scale_low_res = get_best_image_scale(rdmols, opts);
    BOOST_TEST(scale_low_res > 0);
    // lowering the resolution (i.e. a smaller image) should decrease the scale
    BOOST_TEST(scale_low_res < scale_all_three);
    // setting opts.scale shouldn't affect the return value
    opts.scale = 5;
    BOOST_CHECK_CLOSE(get_best_image_scale(rdmols, opts), scale_low_res, 0.01);
}

BOOST_AUTO_TEST_CASE(test_image_scaling)
{
    auto rdmol = rdkit_extensions::to_rdkit("c1nccc2n1ccc2");
    qreal best_scale = get_best_image_scale({rdmol.get()});
    RenderOptions opts;
    opts.scale = best_scale;
    auto bytes_with_best_scale =
        get_image_bytes(*rdmol, ImageFormat::PNG, opts);

    // autoscaling should scale things the same as get_best_image_scale
    auto bytes_with_autoscale = get_image_bytes(*rdmol, ImageFormat::PNG);
    BOOST_TEST(bytes_with_best_scale == bytes_with_autoscale);

    // scale should be ignored if it would make the molecule too big for the
    // image
    opts.scale = best_scale * 2;
    auto bytes_with_scale_too_big =
        get_image_bytes(*rdmol, ImageFormat::PNG, opts);
    BOOST_TEST(bytes_with_best_scale == bytes_with_scale_too_big);

    // for a sanity check, make sure that setting a smaller scale *does* affect
    // the image
    opts.scale = best_scale / 2;
    auto bytes_with_scale_too_small =
        get_image_bytes(*rdmol, ImageFormat::PNG, opts);
    BOOST_TEST(bytes_with_best_scale != bytes_with_scale_too_small);
}
