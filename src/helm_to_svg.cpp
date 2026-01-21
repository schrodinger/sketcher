/* -------------------------------------------------------------------------
 * Simple CLI tool to convert HELM strings to SVG
 * Usage: helm_to_svg <helm_string>
 * Outputs: SVG to stdout
 * ------------------------------------------------------------------------- */

#include <QApplication>
#include <QByteArray>
#include <iostream>
#include <string>

#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/image_generation.h"

using namespace schrodinger;

int main(int argc, char* argv[])
{
    // Qt requires QApplication for rendering
    QApplication app(argc, argv);

    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <helm_string>" << std::endl;
        return 1;
    }

    std::string helm_string = argv[1];

    try {
        // Convert HELM to RDKit molecule
        auto mol = rdkit_extensions::to_rdkit(helm_string,
                                              rdkit_extensions::Format::HELM);

        // Set up rendering options
        sketcher::RenderOptions opts;
        opts.width_height = QSize(400, 400);
        opts.background_color = Qt::white;

        // Generate SVG
        QByteArray svg_bytes =
            sketcher::get_image_bytes(*mol, sketcher::ImageFormat::SVG, opts);

        // Output to stdout
        std::cout << svg_bytes.toStdString();
        return 0;

    } catch (const std::exception& e) {
        std::cerr << "Error generating image for HELM: " << helm_string
                  << ". Error: " << e.what() << std::endl;
        return 1;
    }
}
