#pragma once

#include <filesystem>
#include <fstream>

#include <QtGlobal>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "qapplication_required_fixture.h"

BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::rdkit_extensions::Format)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::AtomQuery)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::AtomTool)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::DrawTool)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::Element)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::EnumerationTool)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::SelectionTool)

/**
 * Gets full path to a file in the testfiles directory in the source
 * @param filename file found in the testfiles folder
 * @return file contents
 */
std::string read_testfile(const std::string& filename)
{
    auto path = std::filesystem::path(std::getenv("SKETCHER_SOURCE_DIR")) /
                "test" / "testfiles" / filename;
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("File not found: " +
                                 std::filesystem::absolute(path).string());
    }

    std::ifstream fh(path.string());
    return std::string(std::istreambuf_iterator<char>(fh),
                       std::istreambuf_iterator<char>());
}

namespace schrodinger
{
namespace sketcher
{

class TestSketcherWidget : public SketcherWidget
{
  public:
    TestSketcherWidget() : SketcherWidget(){};
    using SketcherWidget::addFromString;
    using SketcherWidget::addTextToMolModel;
    using SketcherWidget::copy;
    using SketcherWidget::cut;
    using SketcherWidget::importText;
    using SketcherWidget::m_mol_model;
    using SketcherWidget::m_scene;
    using SketcherWidget::m_sketcher_model;
    using SketcherWidget::m_ui;
    using SketcherWidget::m_undo_stack;
    using SketcherWidget::m_watermark_item;
    using SketcherWidget::paste;

    // using the system clipboard during tests leads to intermittent test
    // failures on buildbot, so we create our own clipboard
    mutable std::string m_clipboard;
    std::string getClipboardContents() const override
    {
        return m_clipboard;
    }
    void setClipboardContents(std::string text) const override
    {
        m_clipboard = text;
    }
};

class TestScene : public Scene
{
  public:
    TestScene(MolModel* mol_model, SketcherModel* sketcher_model,
              QWidget* parent = nullptr) :
        Scene(mol_model, sketcher_model, parent){};
    using Scene::m_atom_to_atom_item;
    using Scene::m_bond_to_bond_item;
    using Scene::m_fonts;
    using Scene::m_mol_model;
    using Scene::m_selection_highlighting_item;
    using Scene::m_sketcher_model;
    using Scene::mousePressEvent;
    using Scene::mouseReleaseEvent;

    static std::shared_ptr<TestScene> getScene()
    {
        auto undo_stack = new QUndoStack();
        auto mol_model = new MolModel(undo_stack);
        auto sketcher_model = new SketcherModel();
        auto scene = std::make_shared<TestScene>(mol_model, sketcher_model);
        undo_stack->setParent(mol_model);
        mol_model->setParent(scene.get());
        sketcher_model->setParent(scene.get());
        return scene;
    }
};

// Helper functions to interfacing with the MolModel through serialized text

void import_mol_text(
    MolModel* mol_model, const std::string& text,
    rdkit_extensions::Format format = rdkit_extensions::Format::AUTO_DETECT)
{
    mol_model->addMol(*to_rdkit(text, format));
}

void import_reaction_text(
    MolModel* mol_model, const std::string& text,
    rdkit_extensions::Format format = rdkit_extensions::Format::AUTO_DETECT)
{
    mol_model->addReaction(*to_rdkit_reaction(text, format));
}

std::string get_mol_text(MolModel* mol_model, rdkit_extensions::Format format)
{
    return rdkit_extensions::to_string(*mol_model->getMolForExport(), format);
}

std::string get_reaction_text(MolModel* mol_model,
                              rdkit_extensions::Format format)
{
    return rdkit_extensions::to_string(*mol_model->getReactionForExport(),
                                       format);
}

} // namespace sketcher
} // namespace schrodinger