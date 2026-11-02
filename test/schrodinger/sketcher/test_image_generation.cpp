/* -------------------------------------------------------------------------
 * Tests schrodinger::sketcher:: rendering APIs
 --------------------------------------------------------------------------- */

#define BOOST_TEST_MODULE sketcher_image_generation

#include <algorithm>
#include <cmath>
#include <vector>

#include <rdkit/GraphMol/GraphMol.h>
#include <QGraphicsItem>
#include <QGraphicsSimpleTextItem>
#include <QList>
#include <QLineF>
#include <QPainter>
#include <QSvgRenderer>
#include <QXmlStreamReader>
#include <QXmlStreamWriter>
#include <boost/filesystem.hpp>
#include <boost/test/framework.hpp>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/helm/monomer_coordgen.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/molviewer/abstract_atom_or_monomer_item.h"
#include "schrodinger/sketcher/molviewer/abstract_bond_or_connector_item.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/monomer_attachment_point_labels.h"
#include "shared_test_helpers.h"
#include "test_common.h"

using boost::unit_test::framework::current_test_case;
using namespace schrodinger::sketcher;
using namespace schrodinger;

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
BOOST_TEST_DONT_PRINT_LOG_VALUE(QSize);

template <typename T> void image_gen_APIs(const T& input)
{
    RenderOptions opts;
    opts.width_height = {400, 400};
    opts.background_color = QColor("#d4e6f1");

    auto qpict = get_qpicture(input, opts);
    BOOST_TEST(qpict.size() > 0);
    BOOST_TEST(qpict.boundingRect().size() == QSize(400, 400));

    auto image = get_qimage(input, opts);
    BOOST_TEST(image.sizeInBytes() > 0);
    BOOST_TEST(image.size() == QSize(400, 400));

    auto bytes = get_image_bytes(input, ImageFormat::PNG, opts);
    BOOST_TEST(bytes.size() > 0);

    const auto tmp_file = "tmp_for_" + current_test_case().p_name.get();
    for (const auto& ext : {".png", ".svg"}) {
        auto outfile = tmp_file + ext;
        save_image_file(input, outfile, opts);
        BOOST_TEST(boost::filesystem::exists(outfile));
    }
    BOOST_REQUIRE_THROW(save_image_file(input, tmp_file + ".invalid", opts),
                        std::out_of_range);
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

/**
 * Make sure that we can generate images containing monomers that don't appear
 * in the monomer database.
 */
BOOST_AUTO_TEST_CASE(test_image_gen_unknown_monomer)
{
    const std::string helm = "PEPTIDE1{[missingMonomer]}$$$$V2.0";
    auto mol = rdkit_extensions::to_rdkit(helm);
    RenderOptions opts;
    BOOST_TEST(!get_qimage(*mol, opts).isNull());
    BOOST_TEST(!get_qimage(helm, opts).isNull());
}

/**
 * Check actual bead bounds and model backbone clearance for a horizontal
 * peptide. Return its edge-to-edge gaps for comparison across label sizes.
 */
static std::vector<qreal> check_peptide_geometry(const TestScene& scene)
{
    const auto mol = scene.m_mol_model->getMol();
    const auto& conformer = mol->getConformer();
    std::vector<QRectF> bounds;
    for (const auto* atom : mol->atoms()) {
        const auto* item = scene.getGraphicsItemForAtom(atom);
        BOOST_CHECK_SMALL(
            QLineF(item->scenePos(),
                   to_scene_xy(conformer.getAtomPos(atom->getIdx())))
                .length(),
            1e-6);
        bounds.push_back(item->sceneBoundingRect());
    }
    for (size_t i = 0; i < bounds.size(); ++i) {
        for (size_t j = i + 1; j < bounds.size(); ++j) {
            BOOST_TEST_CONTEXT("beads " << i << " and " << j)
            {
                BOOST_CHECK(!bounds[i].intersects(bounds[j]));
            }
        }
    }

    std::vector<qreal> gaps;
    const auto connectors =
        scene.getInteractiveItems(InteractiveItemFlag::MONOMER_CONNECTOR);
    BOOST_REQUIRE_EQUAL(connectors.size(), mol->getNumBonds());
    for (const auto* item : connectors) {
        const auto* connector =
            static_cast<const AbstractBondOrConnectorItem*>(item);
        const auto* bond = connector->getBond();
        const auto& start = bounds[bond->getBeginAtomIdx()];
        const auto& end = bounds[bond->getEndAtomIdx()];
        const QLineF line(
            to_scene_xy(conformer.getAtomPos(bond->getBeginAtomIdx())),
            to_scene_xy(conformer.getAtomPos(bond->getEndAtomIdx())));
        BOOST_CHECK_SMALL(line.dy(), 1e-6);
        const auto gap = std::max(start.left(), end.left()) -
                         std::min(start.right(), end.right());
        BOOST_CHECK_GT(gap, std::min(start.height(), end.height()) * 0.25);
        gaps.push_back(gap);

        // The hover tool uses this same factory. Inspect actual text items so
        // sufficient connector length alone cannot hide an AP label collision.
        auto labels = create_attachment_point_labels_for_connector(
            bond, false, BOUND_AP_LABEL_COLOR, scene.m_fonts, &scene);
        BOOST_CHECK_EQUAL(labels.size(), 2);
        for (const auto* label : labels) {
            const auto rect = get_label_ink_rect(
                *static_cast<const QGraphicsSimpleTextItem*>(label));
            for (const auto& bead : bounds) {
                BOOST_CHECK(!rect.intersects(bead));
            }
        }
        if (labels.size() == 2) {
            BOOST_CHECK(
                !get_label_ink_rect(
                     *static_cast<const QGraphicsSimpleTextItem*>(labels[0]))
                     .intersects(get_label_ink_rect(
                         *static_cast<const QGraphicsSimpleTextItem*>(
                             labels[1]))));
        }
        qDeleteAll(labels);
    }
    std::sort(gaps.begin(), gaps.end());
    return gaps;
}

BOOST_AUTO_TEST_CASE(test_initial_large_monomer_labels_expand_layout)
{
    // Inline monomers avoid dependence on an external monomer database.
    const std::string inline_monomer = "[O=C([C@H](CSCCN[H:3])N[H:1])[OH:2]]";
    std::vector<qreal> standard_gaps;
    for (int long_count : {0, 1, 2}) {
        BOOST_TEST_CONTEXT("long label count " << long_count)
        {
            std::string helm = "PEPTIDE1{A.";
            for (int i = 0; i < 2; ++i) {
                helm += (i < long_count ? inline_monomer : "A") + ".";
            }
            helm += "A.A}$$$$V2.0";
            auto scene = TestScene::getScene();
            import_mol_text(scene->m_mol_model, helm);
            BOOST_REQUIRE_EQUAL(scene->m_mol_model->getMol()->getNumAtoms(), 5);
            const auto exported = get_mol_text(scene->m_mol_model,
                                               rdkit_extensions::Format::HELM);
            for (int cleanup_count : {0, 1, 2}) {
                BOOST_TEST_CONTEXT("cleanup count " << cleanup_count)
                {
                    const auto previous = scene->m_mol_model->getMol()
                                              ->getConformer()
                                              .getPositions();
                    if (cleanup_count > 0) {
                        scene->m_mol_model->regenerateCoordinates();
                    }
                    const auto gaps = check_peptide_geometry(*scene);
                    if (long_count == 0 && cleanup_count == 0) {
                        standard_gaps = gaps;
                    }
                    BOOST_REQUIRE_EQUAL(gaps.size(), standard_gaps.size());
                    for (size_t i = 0; i < gaps.size(); ++i) {
                        BOOST_CHECK_GT(gaps[i], standard_gaps[i] * 0.8);
                        BOOST_CHECK_LT(gaps[i], standard_gaps[i] * 1.2);
                    }
                    BOOST_CHECK_EQUAL(
                        get_mol_text(scene->m_mol_model,
                                     rdkit_extensions::Format::HELM),
                        exported);
                    if (cleanup_count == 2) {
                        const auto& current = scene->m_mol_model->getMol()
                                                  ->getConformer()
                                                  .getPositions();
                        for (size_t i = 0; i < current.size(); ++i) {
                            BOOST_CHECK_SMALL(
                                (current[i] - previous[i]).length(), 1e-6);
                        }
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE(test_get_image_bytes_from_scene)
{
    RenderOptions opts;
    auto test_scene = TestScene::getScene();
    import_mol_text(test_scene->m_mol_model, "CCCCC");
    BOOST_TEST(get_image_bytes(*test_scene, ImageFormat::PNG, opts).size() > 0);

    test_scene->m_mol_model->clear();
    import_reaction_text(test_scene->m_mol_model, "CC.CC>>CC");
    BOOST_TEST(get_image_bytes(*test_scene, ImageFormat::PNG, opts).size() > 0);
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
    opts.rdatom_index_to_halo_color = highlights;

    QHash<int, QColor> highlights_2;
    QColor yellow(255, 255, 0);
    highlights[0] = yellow;
    highlights[1] = yellow;

    opts.rdbond_index_to_halo_color = highlights_2;

    auto qimage = get_qimage(*rdmol, opts);
    BOOST_TEST(!qimage.allGray());
}

BOOST_AUTO_TEST_CASE(test_atom_label_halo_color_and_opacity)
{
    for (const auto* smiles : {"O", "N", "S", "[O-]", "[Mg+2]", "[Ne]"}) {
        for (int font_size : {18, 36}) {
            for (bool svg : {false, true}) {
                BOOST_TEST_CONTEXT(smiles << ", font size " << font_size
                                          << ", SVG " << svg)
                {
                    RenderOptions opts;
                    opts.font_size = font_size;
                    opts.scale = 2;
                    auto render = [&](RenderOptions options,
                                      bool hide_labels = false) {
                        if (!svg) {
                            if (hide_labels) {
                                options.rdatom_index_to_line_color[0] =
                                    Qt::transparent;
                            }
                            return get_qimage(smiles, options);
                        }
                        auto svg_data =
                            get_image_bytes(smiles, ImageFormat::SVG, options);
                        if (hide_labels) {
                            // Avoid QRegularExpression here: PCRE2's JIT can
                            // read past the SVG buffer under Valgrind.
                            QXmlStreamReader reader(svg_data);
                            QByteArray without_labels;
                            QXmlStreamWriter writer(&without_labels);
                            while (!reader.atEnd()) {
                                reader.readNext();
                                if (reader.isStartElement() &&
                                    reader.name() == QLatin1String("text")) {
                                    reader.skipCurrentElement();
                                } else {
                                    writer.writeCurrentToken(reader);
                                }
                            }
                            BOOST_REQUIRE_MESSAGE(
                                !reader.hasError(),
                                reader.errorString().toStdString());
                            BOOST_REQUIRE(!writer.hasError());
                            svg_data = without_labels;
                        }
                        QSvgRenderer renderer(svg_data);
                        BOOST_REQUIRE(renderer.isValid());
                        QImage image(options.width_height,
                                     QImage::Format_ARGB32_Premultiplied);
                        image.fill(Qt::transparent);
                        QPainter painter(&image);
                        renderer.render(&painter);
                        return image;
                    };
                    const auto label = render(opts);
                    const auto decorations = render(opts, true);
                    const QColor color(80, 160, 240, 128);
                    opts.rdatom_index_to_halo_color[0] = color;
                    const auto highlighted = render(opts);
                    // Hide the ink to inspect the halo directly beneath every
                    // glyph, including the hydrogen count and charge text.
                    const auto halo = render(opts, true);
                    int label_pixels = 0;
                    int uncovered_pixels = 0;
                    int changed_ink_pixels = 0;
                    for (int y = 0; y < label.height(); ++y) {
                        for (int x = 0; x < label.width(); ++x) {
                            const auto ink = label.pixelColor(x, y);
                            // Valence-error backgrounds are painted above the
                            // halo and are independent of the label's ink.
                            if (ink.alpha() == 0 ||
                                decorations.pixelColor(x, y).alpha() != 0) {
                                continue;
                            }
                            ++label_pixels;
                            const auto background = halo.pixelColor(x, y);
                            // SVG serializes opacity as a decimal, so allow
                            // one unit of rounding in each channel.
                            if (std::abs(background.alpha() - color.alpha()) >
                                    1 ||
                                std::abs(background.red() - color.red()) > 1 ||
                                std::abs(background.green() - color.green()) >
                                    1 ||
                                std::abs(background.blue() - color.blue()) >
                                    1) {
                                ++uncovered_pixels;
                            }
                            if (ink.alpha() == 255 &&
                                highlighted.pixelColor(x, y) != ink) {
                                ++changed_ink_pixels;
                            }
                        }
                    }
                    if (uncovered_pixels > 0) {
                        const auto prefix = QString("tmp_halo_%1_%2_%3_")
                                                .arg(smiles)
                                                .arg(font_size)
                                                .arg(svg);
                        label.save(prefix + "label.png");
                        highlighted.save(prefix + "highlighted.png");
                        halo.save(prefix + "halo.png");
                    }
                    BOOST_CHECK_GT(label_pixels, 0);
                    BOOST_CHECK_EQUAL(uncovered_pixels, 0);
                    BOOST_CHECK_EQUAL(changed_ink_pixels, 0);
                }
            }
        }
    }
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

BOOST_AUTO_TEST_CASE(test_SVG_size_output)
{
    // test that the in SVG output width and height are correctly set in px
    auto rdmol = rdkit_extensions::to_rdkit("c1nccc2n1ccc2");
    RenderOptions opts;
    auto sizes = std::vector<QSize>{{400, 400}, {200, 400}, {300, 200}};
    for (auto size : sizes) {
        opts.width_height = size;
        auto bytes = get_image_bytes(*rdmol, ImageFormat::SVG, opts);
        BOOST_TEST(bytes.size() > 0);
        auto svg = QString::fromUtf8(bytes);
        QString width_height_string =
            QString("svg width=\"%1px\" height=\"%2px\"")
                .arg(size.width())
                .arg(size.height());
        BOOST_TEST(svg.contains(width_height_string));
    }
}

/**
 * Make sure that the trim_image option correctly trims excess white space
 */
BOOST_AUTO_TEST_CASE(test_trim_image)
{
    auto mol = rdkit_extensions::to_rdkit("C1=CC=CC=C1");
    RenderOptions opts;
    opts.width_height = {200, 2000};

    // first, generate the picture/image without trimming and make sure that the
    // images are full size
    auto qpict = get_qpicture(*mol, opts);
    BOOST_TEST(qpict.boundingRect().size() == QSize(200, 2000));
    auto image = get_qimage(*mol, opts);
    BOOST_TEST(image.size() == QSize(200, 2000));

    // enable trimming and make sure that the excess horizontal space is removed
    opts.trim_image = true;
    auto trimmed_wide_qpict = get_qpicture(*mol, opts);
    auto picture_size = trimmed_wide_qpict.boundingRect().size();
    BOOST_TEST(picture_size.width() == 200);
    // the molecule is square-ish, but we don't know (or care) the exact
    // width:height ratio so we just make sure that the height is something
    // vaguely reasonable
    BOOST_TEST(picture_size.height() <= 400);
    BOOST_TEST(picture_size.height() >= 100);

    auto trimmed_wide_image = get_qimage(*mol, opts);
    auto image_size = trimmed_wide_image.size();
    BOOST_TEST(image_size.width() == 200);
    BOOST_TEST(image_size.height() <= 400);
    BOOST_TEST(image_size.height() >= 100);

    // make the requested image size really tall instead of really wide and make
    // sure that the excess vertical space is removed
    opts.width_height = {1000, 100};
    auto trimmed_tall_qpict = get_qpicture(*mol, opts);
    picture_size = trimmed_tall_qpict.boundingRect().size();
    BOOST_TEST(picture_size.width() <= 200);
    BOOST_TEST(picture_size.width() >= 50);
    BOOST_TEST(picture_size.height() == 100);

    auto trimmed_tall_image = get_qimage(*mol, opts);
    image_size = trimmed_tall_image.size();
    BOOST_TEST(image_size.width() <= 200);
    BOOST_TEST(image_size.width() >= 50);
    BOOST_TEST(image_size.height() == 100);
}

/**
 * Verify that Qt's default title and description elements are removed from SVG,
 * and that the font specifications are updated to include fallback fonts
 */
BOOST_AUTO_TEST_CASE(test_SVG_post_processing)
{
    auto rdmol = rdkit_extensions::to_rdkit("c1nccc2n1ccc2");
    RenderOptions opts;
    opts.width_height = {400, 400};

    auto bytes = get_image_bytes(*rdmol, ImageFormat::SVG, opts);
    BOOST_TEST(bytes.size() > 0);
    auto svg = QString::fromUtf8(bytes);

    // Verify Qt's default title and description are not present
    BOOST_TEST(!svg.contains("<title>"));
    BOOST_TEST(!svg.contains("</title>"));
    BOOST_TEST(!svg.contains("<desc>"));
    BOOST_TEST(!svg.contains("</desc>"));
    BOOST_TEST(!svg.contains("Qt SVG Document"));
    BOOST_TEST(!svg.contains("Generated with Qt"));

    // make sure that we specify fallback fonts
    BOOST_TEST(
        svg.contains(R"(font-family="Arimo,Helvetica,Arial,sans-serif")"));
}

/**
 * Test that image generation APIs can handle molecules with more than 200
 * heavy atoms, even though the sketcher widget enforces a 200 atom limit for
 * import.
 */
BOOST_AUTO_TEST_CASE(test_image_generation_with_large_molecules)
{
    // Create a molecule with more than 200 heavy atoms
    // A simple polystyrene-like chain with 101 benzene rings
    std::string large_smiles = "c1ccccc1";
    for (int i = 0; i < 100; ++i) {
        large_smiles += "c1ccccc1";
    }
    auto large_mol = rdkit_extensions::to_rdkit(large_smiles);
    BOOST_REQUIRE(large_mol->getNumAtoms() > 200);
    rdkit_extensions::compute2DCoords(*large_mol);

    // Test that all image generation APIs work with this large molecule
    RenderOptions opts;
    opts.width_height = {800, 800};

    // Test get_qpicture
    auto qpict = get_qpicture(*large_mol, opts);
    BOOST_TEST(qpict.size() > 0);
    BOOST_TEST(qpict.boundingRect().size() == QSize(800, 800));

    // Test get_qimage
    auto image = get_qimage(*large_mol, opts);
    BOOST_TEST(image.sizeInBytes() > 0);
    BOOST_TEST(image.size() == QSize(800, 800));

    // Test get_image_bytes for PNG
    auto png_bytes = get_image_bytes(*large_mol, ImageFormat::PNG, opts);
    BOOST_TEST(png_bytes.size() > 0);

    // Test get_image_bytes for SVG
    auto svg_bytes = get_image_bytes(*large_mol, ImageFormat::SVG, opts);
    BOOST_TEST(svg_bytes.size() > 0);

    // Test save_image_file
    const auto tmp_file =
        "tmp_large_mol_" + current_test_case().p_name.get() + ".png";
    save_image_file(*large_mol, tmp_file, opts);
    BOOST_TEST(boost::filesystem::exists(tmp_file));
}
