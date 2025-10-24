/* -------------------------------------------------------------------------
 * Schrodinger 2D Sketcher Application
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/bind.h>
#endif

#include <QApplication>
#include <QFile>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/sketcher_widget.h"

using schrodinger::rdkit_extensions::Format;
using schrodinger::sketcher::ImageFormat;
using schrodinger::sketcher::SketcherWidget;

// For the WebAssembly build, we need to be able to get the sketcher
// instance we are running from a function/static method. We'll use a
// singleton for that.
SketcherWidget& get_sketcher_instance()
{
    static SketcherWidget instance;
    return instance;
}

void sketcher_import_text(const std::string& text)
{
    auto& sk = get_sketcher_instance();
    sk.addFromString(text);
}

std::string sketcher_export_text(Format format)
{
    auto& sk = get_sketcher_instance();
    return sk.getString(format);
}

std::string sketcher_export_image(ImageFormat format)
{
    auto& sk = get_sketcher_instance();
    return sk.getImageBytes(format).toBase64().toStdString();
}

void sketcher_clear()
{
    auto& sk = get_sketcher_instance();
    sk.clear();
}

bool sketcher_is_empty()
{
    auto& sk = get_sketcher_instance();
    return sk.isEmpty();
}

bool sketcher_has_monomers()
{
    auto& sk = get_sketcher_instance();
    auto mol = sk.getRDKitMolecule();
    return schrodinger::rdkit_extensions::isMonomeric(*mol);
}

void sketcher_changed()
{
#ifdef __EMSCRIPTEN__
    EM_ASM({
        if (Module.sketcher_changed_callback) {
            setTimeout(
                function() {
                    if (Module.sketcher_changed_callback) {
                        Module.sketcher_changed_callback();
                    }
                },
                100);
        }
    });
#endif
}

#ifdef __EMSCRIPTEN__
EMSCRIPTEN_BINDINGS(sketcher)
{
    emscripten::enum_<Format>("Format")
        .value("AUTO_DETECT", Format::AUTO_DETECT)
        .value("RDMOL_BINARY_BASE64", Format::RDMOL_BINARY_BASE64)
        .value("SMILES", Format::SMILES)
        .value("EXTENDED_SMILES", Format::EXTENDED_SMILES)
        .value("SMARTS", Format::SMARTS)
        .value("EXTENDED_SMARTS", Format::EXTENDED_SMARTS)
        .value("MDL_MOLV2000", Format::MDL_MOLV2000)
        .value("MDL_MOLV3000", Format::MDL_MOLV3000)
        .value("MAESTRO", Format::MAESTRO)
        .value("INCHI", Format::INCHI)
        .value("INCHI_KEY", Format::INCHI_KEY)
        .value("PDB", Format::PDB)
        .value("MOL2", Format::MOL2)
        .value("XYZ", Format::XYZ)
        .value("MRV", Format::MRV)
        .value("CDXML", Format::CDXML)
        .value("HELM", Format::HELM)
        .value("FASTA_PEPTIDE", Format::FASTA_PEPTIDE)
        .value("FASTA_DNA", Format::FASTA_DNA)
        .value("FASTA_RNA", Format::FASTA_RNA)
        .value("FASTA", Format::FASTA)
        .value("FMP", Format::FMP)
        .value("CUSTOM_ENTITY", Format::CUSTOM_ENTITY);

    emscripten::enum_<ImageFormat>("ImageFormat")
        .value("PNG", ImageFormat::PNG)
        .value("SVG", ImageFormat::SVG);

    emscripten::function("sketcher_import_text", &sketcher_import_text);
    emscripten::function("sketcher_export_text", &sketcher_export_text);
    emscripten::function("sketcher_export_image", &sketcher_export_image);
    emscripten::function("sketcher_clear", &sketcher_clear);
    emscripten::function("sketcher_is_empty", &sketcher_is_empty);
    emscripten::function("sketcher_has_monomers", &sketcher_has_monomers);
    // see sketcher_changed_callback above
}
#endif

void apply_stylesheet(QApplication& app)
{
    QFile styleFile(":resources/schrodinger_livedesign.qss");
    bool success = styleFile.open(QFile::ReadOnly);
    if (!success) {
        throw std::runtime_error("Could not open style sheet file");
    }
    QString style(styleFile.readAll());
    app.setStyleSheet(style);
}

int main(int argc, char** argv)
{
    QApplication application(argc, argv);
#ifdef SKETCHER_STATIC_DEFINE
    Q_INIT_RESOURCE(sketcher);
#endif

#ifdef __EMSCRIPTEN__
    // Only apply this stylesheet for the WASM build
    apply_stylesheet(application);
    auto& sk = get_sketcher_instance();
    // Qt::WA_AlwaysShowToolTips works around QTBUG-94583
    // (https://bugreports.qt.io/browse/QTBUG-94583), which would otherwise
    // prevent tooltips from showing up. (See SKETCH-2565.)
    sk.setAttribute(Qt::WA_AlwaysShowToolTips);
    QObject::connect(&sk, &SketcherWidget::moleculeChanged, &sketcher_changed);
    QObject::connect(&sk, &SketcherWidget::representationChanged,
                     &sketcher_changed);
#else
    SketcherWidget sk;
#endif

    sk.show();
    return application.exec();
}
