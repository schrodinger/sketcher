#pragma once

#include <fstream>

#include <QtGlobal>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/sketcher_model.h"
#include "schrodinger/test/testfiles.h"
#include "test_markers.h"

BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::AtomQuery)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::AtomTool)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::DrawTool)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::Element)

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
};

} // namespace rdkit_extensions
} // namespace schrodinger