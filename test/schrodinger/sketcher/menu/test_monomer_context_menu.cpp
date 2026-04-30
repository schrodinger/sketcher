#define BOOST_TEST_MODULE Test_Sketcher

#include <algorithm>
#include <unordered_set>

#include <boost/test/unit_test.hpp>

#include <QAction>
#include <QObject>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/menu/monomer_context_menu.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Triggering Delete must emit `deleteRequested` carrying exactly the set
 * of atoms the menu was populated with via `setContextItems`. This is the
 * contract `SketcherWidget` relies on to route the deletion to MolModel.
 */
BOOST_AUTO_TEST_CASE(test_delete_emits_signal_with_context_atoms)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A.G}$$$$V2.0");
    const std::unordered_set<const RDKit::Atom*> input = {
        mol->getAtomWithIdx(0), mol->getAtomWithIdx(1)};

    MonomerContextMenu menu;
    std::unordered_set<const RDKit::Atom*> received;
    int call_count = 0;
    QObject::connect(&menu, &MonomerContextMenu::deleteRequested, &menu,
                     [&](auto atoms) {
                         received = std::move(atoms);
                         ++call_count;
                     });

    menu.setContextItems(input, {}, {}, {}, {});
    auto actions = menu.actions();
    auto it = std::ranges::find_if(
        actions, [](auto* a) { return a->text() == "Delete"; });
    BOOST_REQUIRE(it != actions.end());
    (*it)->trigger();

    BOOST_TEST(call_count == 1);
    BOOST_TEST(received == input);
}

} // namespace sketcher
} // namespace schrodinger
