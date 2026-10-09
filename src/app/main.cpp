/* -------------------------------------------------------------------------
 * Schrodinger Sketcher Application
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#include <emscripten/bind.h>
#include <emscripten/val.h>
#else
#include "crash_handler.h"
#endif

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <QApplication>
#include <QEvent>
#include <QFile>
#include <QIcon>
#include <QPointer>
#include <QString>
#include <QStyleHints>
#include <QTimer>
#include <QWidget>

#include "image_generation_from_js.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "sketcher_instance.h"

using schrodinger::rdkit_extensions::Format;
using schrodinger::sketcher::CarbonLabels;
using schrodinger::sketcher::ColorScheme;
using schrodinger::sketcher::ImageFormat;
using schrodinger::sketcher::RenderOptions;
using schrodinger::sketcher::SketcherWidget;
using schrodinger::sketcher::StereoLabels;

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

#ifdef __EMSCRIPTEN__
std::string get_image_bytes_from_text(const std::string& text,
                                      ImageFormat format)
{
    auto image_bytes = schrodinger::sketcher::get_image_bytes(text, format);
    return image_bytes.toBase64().toStdString();
}

std::string get_image_bytes_from_text(const std::string& text,
                                      ImageFormat format,
                                      const emscripten::val& options)
{
    auto image_bytes = schrodinger::sketcher::get_image_bytes(
        text, format, render_options_from_js(options));
    return image_bytes.toBase64().toStdString();
}
#endif

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

// Retained as a no-op for backwards compatibility with external callers;
// ATOMISTIC_OR_MONOMERIC is now the default interface type (SKETCH-2735).
void sketcher_allow_monomeric(bool /* allow_monomeric */)
{
}

#ifdef __EMSCRIPTEN__
using MonomerDefsInsertionResult =
    std::pair<std::vector<std::string>, std::vector<std::string>>;

emscripten::val
string_vector_to_js_array(const std::vector<std::string>& values)
{
    auto array = emscripten::val::array();
    for (unsigned i = 0; i < values.size(); ++i) {
        array.set(i, values[i]);
    }
    return array;
}

emscripten::val
monomer_defs_insertion_result_to_js(const MonomerDefsInsertionResult& result)
{
    auto object = emscripten::val::object();
    object.set("succeeded", string_vector_to_js_array(result.first));
    object.set("failed", string_vector_to_js_array(result.second));
    return object;
}

emscripten::val sketcher_load_custom_monomers(const std::string& json)
{
    auto& db = schrodinger::rdkit_extensions::MonomerDatabase::instance();
    auto result = db.loadMonomersFromJson(json);
    get_sketcher_instance().setLoadMonomerDatabaseVisible(false);

    return monomer_defs_insertion_result_to_js(result);
}
#endif

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
        .value("FASTA", Format::FASTA);

    emscripten::enum_<ImageFormat>("ImageFormat")
        .value("PNG", ImageFormat::PNG)
        .value("SVG", ImageFormat::SVG);

    emscripten::enum_<StereoLabels>("StereoLabels")
        .value("NONE", StereoLabels::NONE)
        .value("KNOWN", StereoLabels::KNOWN)
        .value("ALL", StereoLabels::ALL);

    emscripten::enum_<CarbonLabels>("CarbonLabels")
        .value("NONE", CarbonLabels::NONE)
        .value("TERMINAL", CarbonLabels::TERMINAL)
        .value("ALL", CarbonLabels::ALL);

    emscripten::enum_<ColorScheme>("ColorScheme")
        .value("DEFAULT", ColorScheme::DEFAULT)
        .value("AVALON", ColorScheme::AVALON)
        .value("CDK", ColorScheme::CDK)
        .value("DARK_MODE", ColorScheme::DARK_MODE)
        .value("BLACK_WHITE", ColorScheme::BLACK_WHITE)
        .value("WHITE_BLACK", ColorScheme::WHITE_BLACK);

    emscripten::function("sketcher_import_text", &sketcher_import_text);
    emscripten::function("sketcher_export_text", &sketcher_export_text);
    emscripten::function("sketcher_export_image", &sketcher_export_image);
    emscripten::function(
        "get_image_bytes",
        emscripten::select_overload<std::string(
            const std::string&, ImageFormat)>(&get_image_bytes_from_text));
    emscripten::function(
        "get_image_bytes",
        emscripten::select_overload<std::string(const std::string&, ImageFormat,
                                                const emscripten::val&)>(
            &get_image_bytes_from_text));
    emscripten::function("sketcher_clear", &sketcher_clear);
    emscripten::function("sketcher_is_empty", &sketcher_is_empty);
    emscripten::function("sketcher_has_monomers", &sketcher_has_monomers);
    emscripten::function("sketcher_allow_monomeric", &sketcher_allow_monomeric);
    emscripten::function("sketcher_load_custom_monomers",
                         &sketcher_load_custom_monomers);
    // see sketcher_changed_callback above
    // the e2e test bindings are registered in playwright_test_bridge.cpp
}
#endif

#ifdef __EMSCRIPTEN__
/**
 * Gives browser focus to the canvas of the given top-level widget, so the
 * browser delivers key events to it again. Qt only does this itself when
 * there's no input context, and the WASM platform plugin always creates one.
 */
void focus_canvas(QWidget& window)
{
    auto canvas = emscripten::val::module_property(
        "specialHTMLTargets")["!qtwindow" + std::to_string(window.winId())];
    if (!canvas.isUndefined()) {
        canvas.call<void>("focus");
    }
}

/**
 * Hands activation back to the window underneath a popup once the popup
 * closes.
 *
 * Qt 6.7's WASM platform plugin activates a popup's window when the popup is
 * shown, but doesn't reactivate anything when it's hidden. Key events then
 * keep going to the hidden popup (and, for popups that close themselves, the
 * browser leaves focus on the page body), so shortcuts stop working until the
 * user clicks in the sketcher (SKETCH-1652). This is fixed upstream by
 * QTBUG-144692, which isn't in a Qt release yet; remove this once we're on one
 * that includes it.
 */
class PopupActivationRestorer : public QObject
{
  public:
    using QObject::QObject;

    bool eventFilter(QObject* watched, QEvent* event) override
    {
        auto* popup = qobject_cast<QWidget*>(watched);
        if (event->type() == QEvent::Hide && popup != nullptr &&
            popup->windowType() == Qt::Popup) {
            auto* parent = popup->parentWidget();
            QPointer<QWidget> window =
                parent != nullptr ? parent->window() : &get_sketcher_instance();
            // Defer until Qt has finished closing the popup, and leave things
            // alone if another popup (e.g. a parent menu) is still open or if
            // the popup's action opened a modal dialog
            QTimer::singleShot(0, this, [window]() {
                auto* modal = QApplication::activeModalWidget();
                if (window != nullptr && window->isVisible() &&
                    QApplication::activePopupWidget() == nullptr &&
                    (modal == nullptr || modal == window)) {
                    window->activateWindow();
                    focus_canvas(*window);
                }
            });
        }
        return QObject::eventFilter(watched, event);
    }
};
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
#ifndef __EMSCRIPTEN__
    schrodinger::install_crash_handlers();
#endif

    QApplication application(argc, argv);
#ifdef SKETCHER_STATIC_DEFINE
    Q_INIT_RESOURCE(sketcher);
#endif
    QApplication::setWindowIcon(QIcon(":icons/sketcher-logo.svg"));

    // In Qt 6.8 and newer, Qt will try to automatically apply a dark mode color
    // scheme if the system and/or web browser is set to dark mode. The result
    // looks terrible, so switch back to light mode.
#if QT_VERSION >= QT_VERSION_CHECK(6, 8, 0)
    QApplication::styleHints()->setColorScheme(Qt::ColorScheme::Light);
#endif

#ifdef __EMSCRIPTEN__
    // Only apply this stylesheet for the WASM build
    apply_stylesheet(application);
    auto& sk = get_sketcher_instance();
    // Qt::WA_AlwaysShowToolTips works around QTBUG-94583
    // (https://bugreports.qt.io/browse/QTBUG-94583), which would otherwise
    // prevent tooltips from showing up. (See SKETCH-2565.)
    sk.setAttribute(Qt::WA_AlwaysShowToolTips);
    application.installEventFilter(new PopupActivationRestorer(&application));
    QObject::connect(&sk, &SketcherWidget::moleculeChanged, &sketcher_changed);
    QObject::connect(&sk, &SketcherWidget::representationChanged,
                     &sketcher_changed);
#else
    SketcherWidget sk;
#endif

    sk.show();
    return application.exec();
}
