#pragma once

#include <fstream>

#include <QtGlobal>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/scene.h"
#include "schrodinger/sketcher/rdkit/molops.h"
#include "schrodinger/test/testfiles.h"
#include "test_markers.h"

BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::rdkit_extensions::Format)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::AtomQuery)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::AtomTool)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::DrawTool)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::Element)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::EnumerationTool)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::SelectionTool)

/* Use this fixture to construct global fixture which would provide
 * QApplication for a whole suite of a given test file.
 * It should always be used in the BOOST_GLOBAL_FIXTURE
 * This fixture is not suitable as a global fixture if test depends on
 * QApplication.
 */
class Test_Sketcher_global_fixture : public Test_QAPP_DisplayRequiredFixture
{
  public:
    Test_Sketcher_global_fixture() : Test_QAPP_DisplayRequiredFixture(true)
    {
    }
};

std::string read_testfile(const std::string& filename)
{
    std::ifstream fh(schrodinger::test::mmshare_testfile(filename));
    return std::string(std::istreambuf_iterator<char>(fh),
                       std::istreambuf_iterator<char>());
}

namespace schrodinger
{
namespace rdkit_extensions
{

const std::vector<Format> TEXT_FORMATS = {
    Format::SMILES,       Format::EXTENDED_SMILES,
    Format::SMARTS,       Format::MAESTRO,
    Format::MDL_MOLV2000, Format::MDL_MOLV3000,
    Format::INCHI,        Format::PDB,
    Format::XYZ,          Format::RDMOL_BINARY_BASE64,
};

const std::vector<Format> REACTION_TEXT_FORMATS = {
    Format::SMILES,
    Format::SMARTS,
    Format::MDL_MOLV2000,
    Format::MDL_MOLV3000,
    Format::RDMOL_BINARY_BASE64,
};

} // namespace rdkit_extensions
} // namespace schrodinger

namespace schrodinger
{
namespace sketcher
{

class TestScene : public Scene
{
  public:
    TestScene() : Scene(new MolModel(new QUndoStack()), new SketcherModel()){};
    using Scene::m_atom_item_settings;
    using Scene::m_atom_to_atom_item;
    using Scene::m_bond_item_settings;
    using Scene::m_bond_to_bond_item;
    using Scene::m_fonts;
    using Scene::m_mol_model;
    using Scene::m_selection_highlighting_item;
    using Scene::m_sketcher_model;
};

// Helper functions to interfacing with the MolModel through serialized text

void import_mol_text(
    MolModel* mol_model, const std::string& text,
    rdkit_extensions::Format format = rdkit_extensions::Format::AUTO_DETECT)
{
    mol_model->addMol(*text_to_mol(text, format));
}

void import_reaction_text(
    MolModel* mol_model, const std::string& text,
    rdkit_extensions::Format format = rdkit_extensions::Format::AUTO_DETECT)
{
    mol_model->addReaction(*to_rdkit_reaction(text, format));
}

std::string get_mol_text(MolModel* mol_model, rdkit_extensions::Format format)
{
    return rdkit_extensions::to_string(*mol_model->getMol(), format);
}

std::string get_reaction_text(MolModel* mol_model,
                              rdkit_extensions::Format format)
{
    return rdkit_extensions::to_string(*mol_model->getReaction(), format);
}

} // namespace sketcher
} // namespace schrodinger