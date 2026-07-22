#define BOOST_TEST_MODULE Test_Sketcher

#include <algorithm>
#include <string>
#include <string_view>
#include <unordered_set>
#include <utility>
#include <vector>

#include <boost/test/unit_test.hpp>

#include <QAction>
#include <QMenu>
#include <QObject>
#include <QString>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/ROMol.h>

#include "../../rdkit_extensions/test_common.h"
#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/sketcher/menu/monomer_context_menu.h"
#include "schrodinger/sketcher/rdkit/monomeric.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

class TestMonomerContextMenu : public MonomerContextMenu
{
  public:
    TestMonomerContextMenu() : MonomerContextMenu()
    {
    }
    using MonomerContextMenu::updateActions;
};

namespace
{

QAction* find_action(QMenu& menu, const QString& text)
{
    for (auto* a : menu.actions()) {
        if (a->text() == text) {
            return a;
        }
    }
    return nullptr;
}

std::unordered_set<const RDKit::Atom*>
atoms_from(const RDKit::ROMol& mol, std::initializer_list<unsigned int> idxs)
{
    std::unordered_set<const RDKit::Atom*> out;
    for (auto i : idxs) {
        out.insert(mol.getAtomWithIdx(i));
    }
    return out;
}

// Captures the most recent `mutateMonomerRequested` emission.
struct MutateRequestCapture {
    std::vector<MonomerMutation> last_mutations;
    QString last_description;
    int call_count = 0;

    void connectTo(MonomerContextMenu& menu)
    {
        QObject::connect(&menu, &MonomerContextMenu::mutateMonomerRequested,
                         &menu,
                         [this](auto mutations, const QString& description) {
                             last_mutations = std::move(mutations);
                             last_description = description;
                             ++call_count;
                         });
    }
};

// Captures the most recent `addComplementaryStrandRequested` emission.
struct ComplementRequestCapture {
    std::unordered_set<const RDKit::Atom*> last_bases;
    int call_count = 0;

    void connectTo(MonomerContextMenu& menu)
    {
        QObject::connect(&menu,
                         &MonomerContextMenu::addComplementaryStrandRequested,
                         &menu, [this](auto bases) {
                             last_bases = std::move(bases);
                             ++call_count;
                         });
    }
};

// RAII helper for tests that need to inject custom monomer definitions.
// Pairs with the LocalMonomerDbFixture global fixture, which redirects
// the singleton to a per-binary scratch DB so we never touch the user's
// real custom DB.
struct CustomDbScope {
    explicit CustomDbScope(std::string_view json)
    {
        rdkit_extensions::MonomerDatabase::instance().loadMonomersFromJson(
            json);
    }
    ~CustomDbScope()
    {
        rdkit_extensions::MonomerDatabase::instance().resetMonomerDefinitions();
    }
};

} // namespace

/**
 * Triggering Delete must emit `deleteRequested` carrying exactly the
 * set of atoms the menu was populated with via `setContextItems`.
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
    auto* del = find_action(menu, "Delete");
    BOOST_REQUIRE(del != nullptr);
    del->trigger();

    BOOST_TEST(call_count == 1);
    BOOST_TEST(received == input);
}

/**
 * On a peptide selection, all AA-specific actions are visible alongside
 * Delete. Protonate is permanently disabled. Mutate Residue is a
 * populated submenu with one entry per natural amino acid.
 */
BOOST_AUTO_TEST_CASE(test_actions_for_peptide_selection)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* mutate = find_action(menu, "Mutate Residue");
    auto* d_form = find_action(menu, "Set D-Form");
    auto* protonate = find_action(menu, "Protonate");
    auto* del = find_action(menu, "Delete");

    BOOST_REQUIRE(mutate != nullptr);
    BOOST_REQUIRE(d_form != nullptr);
    BOOST_REQUIRE(protonate != nullptr);
    BOOST_REQUIRE(del != nullptr);

    BOOST_TEST(mutate->isVisible());
    BOOST_TEST(d_form->isVisible());
    BOOST_TEST(protonate->isVisible());
    BOOST_TEST(del->isVisible());
    BOOST_TEST(!protonate->isEnabled());

    BOOST_REQUIRE(mutate->menu() != nullptr);
    BOOST_TEST(mutate->menu()->actions().size() == 21u);
}

/**
 * Non-peptide selections (NA bases / sugars / phosphates / CHEM) only
 * see Delete; AA-specific actions hide themselves.
 */
BOOST_AUTO_TEST_CASE(test_non_peptide_hides_aa_actions)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0, 1, 2}), {}, {}, {}, {});
    menu.updateActions();

    BOOST_TEST(!find_action(menu, "Mutate Residue")->isVisible());
    BOOST_TEST(!find_action(menu, "Protonate")->isVisible());
    // updateActions doesn't touch the label text when !peptide, so it
    // still reads its constructor default "Set D-Form".
    auto* d_form = find_action(menu, "Set D-Form");
    BOOST_REQUIRE(d_form != nullptr);
    BOOST_TEST(!d_form->isVisible());
    BOOST_TEST(find_action(menu, "Delete")->isVisible());
}

/**
 * Empty selection should never enter the AA branch.
 */
BOOST_AUTO_TEST_CASE(test_empty_selection_hides_aa_actions)
{
    TestMonomerContextMenu menu;
    menu.setContextItems({}, {}, {}, {}, {});
    menu.updateActions();

    BOOST_TEST(!find_action(menu, "Mutate Residue")->isVisible());
    BOOST_TEST(find_action(menu, "Delete")->isVisible());
}

/**
 * A Mutate Residue leaf emits a single-pair batch with an empty
 * description (so the receiver falls through to MolModel's default
 * undo label).
 */
BOOST_AUTO_TEST_CASE(test_mutate_residue_leaf_emits_single_pair_batch)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* mutate = find_action(menu, "Mutate Residue");
    BOOST_REQUIRE(mutate->menu() != nullptr);
    auto* asp = find_action(*mutate->menu(), "Aspartate (D)");
    BOOST_REQUIRE(asp != nullptr);
    BOOST_REQUIRE(asp->menu() != nullptr);
    auto* d_leaf = find_action(*asp->menu(), "D");
    BOOST_REQUIRE(d_leaf != nullptr);
    d_leaf->trigger();

    BOOST_TEST(capture.call_count == 1);
    BOOST_TEST(capture.last_description.isEmpty());
    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "D");
    BOOST_REQUIRE(capture.last_mutations[0].atoms.size() == 1u);
    BOOST_TEST(*capture.last_mutations[0].atoms.begin() ==
               mol->getAtomWithIdx(0));
}

/**
 * The Mutate Residue submenu populates non-natural analogs from the DB.
 * "meD" (N-Methyl-L-Aspartic acid) should appear under Aspartate (D) —
 * a sentinel guarding against DB query / sort regressions.
 */
BOOST_AUTO_TEST_CASE(test_mutate_residue_submenu_contains_known_analog)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* mutate = find_action(menu, "Mutate Residue");
    BOOST_REQUIRE(mutate->menu() != nullptr);
    auto* asp = find_action(*mutate->menu(), "Aspartate (D)");
    BOOST_REQUIRE(asp != nullptr);
    BOOST_REQUIRE(asp->menu() != nullptr);

    bool found_meD = false;
    for (auto* a : asp->menu()->actions()) {
        if (a->text().startsWith("meD")) {
            found_meD = true;
            break;
        }
    }
    BOOST_TEST(found_meD);
}

/**
 * The natural amino acid leaf must appear exactly once under its own
 * sub-submenu — we add it explicitly while the DB query
 * (getMonomersByNaturalAnalog) excludes it from the analog list. If
 * that DB-side contract ever changes we'd start seeing a duplicate
 * "D" leaf under "Aspartate (D)"; this test guards against that.
 */
BOOST_AUTO_TEST_CASE(test_natural_aa_appears_exactly_once_in_submenu)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* mutate = find_action(menu, "Mutate Residue");
    BOOST_REQUIRE(mutate->menu() != nullptr);
    auto* asp = find_action(*mutate->menu(), "Aspartate (D)");
    BOOST_REQUIRE(asp != nullptr);
    BOOST_REQUIRE(asp->menu() != nullptr);

    int d_count = 0;
    for (auto* a : asp->menu()->actions()) {
        if (a->text() == "D") {
            ++d_count;
        }
    }
    BOOST_TEST(d_count == 1);
}

/**
 * Mutate Residue → X on a residue already at X must NOT emit, so we
 * don't push a no-op mutation onto the undo stack.
 */
BOOST_AUTO_TEST_CASE(test_mutate_residue_to_current_symbol_is_noop)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{D}$$$$V2.0");
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* mutate = find_action(menu, "Mutate Residue");
    BOOST_REQUIRE(mutate->menu() != nullptr);
    auto* asp = find_action(*mutate->menu(), "Aspartate (D)");
    BOOST_REQUIRE(asp != nullptr);
    auto* d_leaf = find_action(*asp->menu(), "D");
    BOOST_REQUIRE(d_leaf != nullptr);
    d_leaf->trigger();

    BOOST_TEST(capture.call_count == 0);
}

/**
 * Triggering "Set L-Form" on a D residue emits a batch with the
 * stripped (L) symbol and the matching macro label — symmetric with
 * test_d_form_toggle_emits_grouped_batch but covering the D→L path.
 */
BOOST_AUTO_TEST_CASE(test_set_l_form_toggle_emits_l_symbol)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{[dD]}$$$$V2.0");
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* l_form = find_action(menu, "Set L-Form");
    BOOST_REQUIRE(l_form != nullptr);
    BOOST_TEST(l_form->isEnabled());
    l_form->trigger();

    BOOST_TEST(capture.call_count == 1);
    BOOST_TEST(capture.last_description == "Set L-Form");
    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "D");
}

/**
 * The D-/L-form action label flips based on the residue's current form.
 */
BOOST_AUTO_TEST_CASE(test_d_form_label_flips_with_selection)
{
    auto l_mol = rdkit_extensions::to_rdkit("PEPTIDE1{D}$$$$V2.0");
    {
        TestMonomerContextMenu menu;
        menu.setContextItems(atoms_from(*l_mol, {0}), {}, {}, {}, {});
        menu.updateActions();
        BOOST_TEST(find_action(menu, "Set D-Form") != nullptr);
        BOOST_TEST(find_action(menu, "Set L-Form") == nullptr);
    }

    auto d_mol = rdkit_extensions::to_rdkit("PEPTIDE1{[dD]}$$$$V2.0");
    {
        TestMonomerContextMenu menu;
        menu.setContextItems(atoms_from(*d_mol, {0}), {}, {}, {}, {});
        menu.updateActions();
        BOOST_TEST(find_action(menu, "Set L-Form") != nullptr);
        BOOST_TEST(find_action(menu, "Set D-Form") == nullptr);
    }
}

/**
 * Mixed L+D selection defaults to "Set D-Form" — clicking it normalises
 * to D, leaving the existing D residue alone.
 */
BOOST_AUTO_TEST_CASE(test_d_form_mixed_selection_defaults_to_set_d)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{D.[dD]}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0, 1}), {}, {}, {}, {});
    menu.updateActions();

    BOOST_TEST(find_action(menu, "Set D-Form") != nullptr);
    BOOST_TEST(find_action(menu, "Set L-Form") == nullptr);
}

/**
 * Glycine has no D-form counterpart in the DB (it's achiral), so the
 * action must disable rather than silently do nothing on click.
 */
BOOST_AUTO_TEST_CASE(test_d_form_disabled_when_no_db_counterpart)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{G}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* d_form = find_action(menu, "Set D-Form");
    BOOST_REQUIRE(d_form != nullptr);
    BOOST_TEST(!d_form->isEnabled());
}

/**
 * D-Form toggle on a mixed L+D selection emits one batch entry per
 * distinct target symbol (here only the L residue needs mutating).
 * Description carries the macro label that the receiver feeds to
 * QUndoStack::beginMacro.
 */
BOOST_AUTO_TEST_CASE(test_d_form_toggle_emits_grouped_batch)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{D.[dD]}$$$$V2.0");
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0, 1}), {}, {}, {}, {});
    menu.updateActions();

    auto* d_form = find_action(menu, "Set D-Form");
    BOOST_REQUIRE(d_form != nullptr);
    BOOST_TEST(d_form->isEnabled());
    d_form->trigger();

    BOOST_TEST(capture.call_count == 1);
    BOOST_TEST(capture.last_description == "Set D-Form");
    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "dD");
    BOOST_REQUIRE(capture.last_mutations[0].atoms.size() == 1u);
    BOOST_TEST(*capture.last_mutations[0].atoms.begin() ==
               mol->getAtomWithIdx(0));
}

/**
 * A custom monomer named "dXyz" whose NATURAL_ANALOG points at an
 * unrelated residue ("K") must NOT be classified as the D-form of
 * "Xyz" — the prefix is meaningless without DB corroboration. The
 * Set D-Form action stays disabled (no toggle target) and labels as
 * "Set D-Form" (not "Set L-Form").
 */
BOOST_AUTO_TEST_CASE(test_custom_d_prefix_with_mismatched_analog_is_not_d_form)
{
    constexpr std::string_view json = R"([{
        "symbol": "dXyz",
        "polymer_type": "PEPTIDE",
        "natural_analog": "K",
        "smiles": "C[C@@H](N[H:1])C(=O)[OH:2]",
        "core_smiles": "C[C@@H](N[H:1])C(=O)[OH:2]",
        "name": "Custom dXyz",
        "monomer_type": "Backbone",
        "author": "test",
        "pdbcode": "DXYZ"
    }])";
    CustomDbScope scope(json);

    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{[dXyz]}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* d_form = find_action(menu, "Set D-Form");
    BOOST_REQUIRE(d_form != nullptr);
    BOOST_TEST(!d_form->isEnabled());
    BOOST_TEST(find_action(menu, "Set L-Form") == nullptr);
}

/**
 * A realistic custom L+D pair: alpha-aminobutyric acid ("Abu") and its
 * D-form ("dAbu"), both tagged as Ala-class variants
 * (NATURAL_ANALOG="A"). Two symbols with the same NATURAL_ANALOG and
 * the `d` prefix relationship are recognised as a stereo pair, so the
 * toggle enables and emits "dAbu" on click.
 */
BOOST_AUTO_TEST_CASE(test_custom_well_formed_d_form_pair_enables_toggle)
{
    constexpr std::string_view json = R"([
        {
            "symbol": "Abu",
            "polymer_type": "PEPTIDE",
            "natural_analog": "A",
            "smiles": "CC[C@H](N[H:1])C(=O)[OH:2]",
            "core_smiles": "CC[C@H](N[H:1])C(=O)[OH:2]",
            "name": "alpha-Aminobutyric acid",
            "monomer_type": "Backbone",
            "author": "test",
            "pdbcode": "ABU"
        },
        {
            "symbol": "dAbu",
            "polymer_type": "PEPTIDE",
            "natural_analog": "A",
            "smiles": "CC[C@@H](N[H:1])C(=O)[OH:2]",
            "core_smiles": "CC[C@@H](N[H:1])C(=O)[OH:2]",
            "name": "D-alpha-Aminobutyric acid",
            "monomer_type": "Backbone",
            "author": "test",
            "pdbcode": "DAB"
        }
    ])";
    CustomDbScope scope(json);

    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{[Abu]}$$$$V2.0");
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* d_form = find_action(menu, "Set D-Form");
    BOOST_REQUIRE(d_form != nullptr);
    BOOST_TEST(d_form->isEnabled());
    d_form->trigger();
    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "dAbu");
}

/**
 * Both symbols exist and the `d` prefix relationship looks right, but
 * their NATURAL_ANALOGs disagree — the author has tagged them as
 * variants of different natural amino acids, so they aren't actually
 * stereo partners. The classifier rejects the pair.
 */
BOOST_AUTO_TEST_CASE(test_custom_paired_with_different_analogs_is_not_d_form)
{
    constexpr std::string_view json = R"([
        {
            "symbol": "Foo",
            "polymer_type": "PEPTIDE",
            "natural_analog": "A",
            "smiles": "C[C@H](N[H:1])C(=O)[OH:2]",
            "core_smiles": "C[C@H](N[H:1])C(=O)[OH:2]",
            "name": "Custom Foo",
            "monomer_type": "Backbone",
            "author": "test",
            "pdbcode": "FOO"
        },
        {
            "symbol": "dFoo",
            "polymer_type": "PEPTIDE",
            "natural_analog": "K",
            "smiles": "C[C@@H](N[H:1])C(=O)[OH:2]",
            "core_smiles": "C[C@@H](N[H:1])C(=O)[OH:2]",
            "name": "Custom dFoo",
            "monomer_type": "Backbone",
            "author": "test",
            "pdbcode": "DFO"
        }
    ])";
    CustomDbScope scope(json);

    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{[Foo]}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* d_form = find_action(menu, "Set D-Form");
    BOOST_REQUIRE(d_form != nullptr);
    BOOST_TEST(!d_form->isEnabled());
}

/**
 * Base selection shows Mutate Base; AA and sugar actions hide. Sanity
 * check on the dispatcher in updateActions().
 */
BOOST_AUTO_TEST_CASE(test_na_base_selection_shows_mutate_base_only)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    // R = atom 0, A (base) = atom 1, P = atom 2 in this molecule.
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {1}), {}, {}, {}, {},
                         mol->getAtomWithIdx(1));
    menu.updateActions();

    BOOST_TEST(find_action(menu, "Mutate Base")->isVisible());
    BOOST_TEST(!find_action(menu, "Mutate Residue")->isVisible());
    BOOST_TEST(!find_action(menu, "Set D-Form")->isVisible());
    auto* sugar = find_action(menu, "Change R to dR");
    BOOST_REQUIRE(sugar != nullptr);
    BOOST_TEST(!sugar->isVisible());
    BOOST_TEST(find_action(menu, "Delete")->isVisible());
}

/**
 * Sugar (R) selection: the toggle reads "Change R to dR" and is enabled
 * when the right-clicked atom is the only sugar in the selection.
 */
BOOST_AUTO_TEST_CASE(test_sugar_toggle_label_R_to_dR)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {},
                         mol->getAtomWithIdx(0));
    menu.updateActions();

    auto* toggle = find_action(menu, "Change R to dR");
    BOOST_REQUIRE(toggle != nullptr);
    BOOST_TEST(toggle->isVisible());
    BOOST_TEST(toggle->isEnabled());
    BOOST_TEST(find_action(menu, "Mutate Base") != nullptr);
    BOOST_TEST(!find_action(menu, "Mutate Base")->isVisible());
}

/**
 * Sugar (dR) selection: the toggle reads "Change dR to R".
 */
BOOST_AUTO_TEST_CASE(test_sugar_toggle_label_dR_to_R)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{[dR](T)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {},
                         mol->getAtomWithIdx(0));
    menu.updateActions();

    auto* toggle = find_action(menu, "Change dR to R");
    BOOST_REQUIRE(toggle != nullptr);
    BOOST_TEST(toggle->isVisible());
    BOOST_TEST(toggle->isEnabled());
}

/**
 * Mixed {R, dR} selection: direction comes from m_primary_atom. With
 * primary = R, all sugars become dR. The atom already at dR is excluded
 * from the emitted batch (no no-op).
 */
BOOST_AUTO_TEST_CASE(test_sugar_toggle_mixed_selection_primary_R)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P.[dR](T)P}$$$$V2.0");
    // R = atom 0, A = atom 1, P = atom 2, dR = atom 3, T = atom 4, P = atom 5.
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0, 3}), {}, {}, {}, {},
                         mol->getAtomWithIdx(0));
    menu.updateActions();

    auto* toggle = find_action(menu, "Change R to dR");
    BOOST_REQUIRE(toggle != nullptr);
    BOOST_TEST(toggle->isEnabled());
    toggle->trigger();

    BOOST_TEST(capture.call_count == 1);
    BOOST_TEST(capture.last_description == "Change R to dR");
    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "dR");
    BOOST_REQUIRE(capture.last_mutations[0].atoms.size() == 1u);
    BOOST_TEST(*capture.last_mutations[0].atoms.begin() ==
               mol->getAtomWithIdx(0));
}

/**
 * Same selection but primary = dR: direction reverses, all sugars become R.
 */
BOOST_AUTO_TEST_CASE(test_sugar_toggle_mixed_selection_primary_dR)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P.[dR](T)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0, 3}), {}, {}, {}, {},
                         mol->getAtomWithIdx(3));
    menu.updateActions();

    auto* toggle = find_action(menu, "Change dR to R");
    BOOST_REQUIRE(toggle != nullptr);
    BOOST_TEST(toggle->isEnabled());
    toggle->trigger();

    BOOST_TEST(capture.call_count == 1);
    BOOST_TEST(capture.last_description == "Change dR to R");
    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "R");
    BOOST_REQUIRE(capture.last_mutations[0].atoms.size() == 1u);
    BOOST_TEST(*capture.last_mutations[0].atoms.begin() ==
               mol->getAtomWithIdx(3));
}

/**
 * No primary atom (e.g. menu shown via API outside the right-click flow):
 * direction is ambiguous, action disabled.
 */
BOOST_AUTO_TEST_CASE(test_sugar_toggle_disabled_when_no_primary_atom)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* toggle = find_action(menu, "Change R to dR");
    BOOST_REQUIRE(toggle != nullptr);
    BOOST_TEST(toggle->isVisible());
    BOOST_TEST(!toggle->isEnabled());
}

/**
 * The Mutate Base submenu has one top-level entry per natural base
 * (A/C/G/U/T) in display order; each contains at least the natural
 * symbol leaf.
 */
BOOST_AUTO_TEST_CASE(test_mutate_base_submenu_structure)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {1}), {}, {}, {}, {},
                         mol->getAtomWithIdx(1));
    menu.updateActions();

    auto* mutate = find_action(menu, "Mutate Base");
    BOOST_REQUIRE(mutate != nullptr);
    BOOST_REQUIRE(mutate->menu() != nullptr);
    BOOST_TEST(mutate->menu()->actions().size() == 5u);

    for (auto* expected : {"Adenine (A)", "Cytosine (C)", "Guanine (G)",
                           "Uracil (U)", "Thymine (T)"}) {
        auto* sub = find_action(*mutate->menu(), expected);
        BOOST_REQUIRE(sub != nullptr);
        BOOST_REQUIRE(sub->menu() != nullptr);
        BOOST_TEST(sub->menu()->actions().size() >= 1u);
    }
}

/**
 * Picking a different base emits a single-pair batch with empty
 * description (per the same convention as Mutate Residue leaves).
 */
BOOST_AUTO_TEST_CASE(test_mutate_base_emits_single_pair_batch)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {1}), {}, {}, {}, {},
                         mol->getAtomWithIdx(1));
    menu.updateActions();

    auto* mutate = find_action(menu, "Mutate Base");
    auto* thy = find_action(*mutate->menu(), "Thymine (T)");
    BOOST_REQUIRE(thy != nullptr);
    auto* t_leaf = find_action(*thy->menu(), "T");
    BOOST_REQUIRE(t_leaf != nullptr);
    t_leaf->trigger();

    BOOST_TEST(capture.call_count == 1);
    BOOST_TEST(capture.last_description.isEmpty());
    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "T");
    BOOST_REQUIRE(capture.last_mutations[0].atoms.size() == 1u);
    BOOST_TEST(*capture.last_mutations[0].atoms.begin() ==
               mol->getAtomWithIdx(1));
}

/**
 * Mutate Base → X on a base already at X must NOT emit, mirroring the
 * Mutate Residue no-op guard.
 */
BOOST_AUTO_TEST_CASE(test_mutate_base_to_current_symbol_is_noop)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {1}), {}, {}, {}, {},
                         mol->getAtomWithIdx(1));
    menu.updateActions();

    auto* mutate = find_action(menu, "Mutate Base");
    auto* ade = find_action(*mutate->menu(), "Adenine (A)");
    BOOST_REQUIRE(ade != nullptr);
    auto* a_leaf = find_action(*ade->menu(), "A");
    BOOST_REQUIRE(a_leaf != nullptr);
    a_leaf->trigger();

    BOOST_TEST(capture.call_count == 0);
}

/**
 * Custom sugar with mismatched NATURAL_ANALOG must not be classified as
 * the deoxy form — toggle stays disabled. Parallel of the AA case.
 */
BOOST_AUTO_TEST_CASE(
    test_custom_d_prefix_sugar_with_mismatched_analog_is_not_deoxy)
{
    constexpr std::string_view json = R"([{
        "symbol": "dXur",
        "polymer_type": "RNA",
        "natural_analog": "K",
        "smiles": "OC[C@H]1O[C@H]([H:1])[C@@H]([OH:2])[C@@H]1O",
        "core_smiles": "OC[C@H]1O[C@H]([H:1])[C@@H]([OH:2])[C@@H]1O",
        "name": "Custom dXur",
        "monomer_type": "Backbone",
        "author": "test",
        "pdbcode": "DXR"
    }])";
    CustomDbScope scope(json);

    auto mol = rdkit_extensions::to_rdkit("RNA1{[dXur]}$$$$V2.0");
    // Classifier precondition — name must end in 'r' to be NA_SUGAR.
    BOOST_REQUIRE(get_monomer_type(mol->getAtomWithIdx(0)) ==
                  MonomerType::NA_SUGAR);

    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {},
                         mol->getAtomWithIdx(0));
    menu.updateActions();

    // No DB counterpart → falls back to default label, disabled.
    auto* toggle = find_action(menu, "Change R to dR");
    BOOST_REQUIRE(toggle != nullptr);
    BOOST_TEST(!toggle->isEnabled());
}

/**
 * Custom sugar pair "Xur"/"dXur" with matching NATURAL_ANALOG: toggle
 * enables and emits the deoxy symbol on click.
 */
BOOST_AUTO_TEST_CASE(test_custom_well_formed_sugar_pair_enables_toggle)
{
    constexpr std::string_view json = R"([
        {
            "symbol": "Xur",
            "polymer_type": "RNA",
            "natural_analog": "Xur",
            "smiles": "OC[C@H]1O[C@H]([H:1])[C@@H]([OH:2])[C@@H]1O",
            "core_smiles": "OC[C@H]1O[C@H]([H:1])[C@@H]([OH:2])[C@@H]1O",
            "name": "Custom Xur",
            "monomer_type": "Backbone",
            "author": "test",
            "pdbcode": "XUR"
        },
        {
            "symbol": "dXur",
            "polymer_type": "RNA",
            "natural_analog": "Xur",
            "smiles": "OC[C@H]1O[C@H]([H:1])C[C@@H]1[OH:2]",
            "core_smiles": "OC[C@H]1O[C@H]([H:1])C[C@@H]1[OH:2]",
            "name": "Custom dXur",
            "monomer_type": "Backbone",
            "author": "test",
            "pdbcode": "DXR"
        }
    ])";
    CustomDbScope scope(json);

    auto mol = rdkit_extensions::to_rdkit("RNA1{[Xur]}$$$$V2.0");
    // Classifier precondition — name must end in 'r' to be NA_SUGAR.
    BOOST_REQUIRE(get_monomer_type(mol->getAtomWithIdx(0)) ==
                  MonomerType::NA_SUGAR);

    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {},
                         mol->getAtomWithIdx(0));
    menu.updateActions();

    auto* toggle = find_action(menu, "Change Xur to dXur");
    BOOST_REQUIRE(toggle != nullptr);
    BOOST_TEST(toggle->isEnabled());
    toggle->trigger();

    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "dXur");
}

/**
 * Sugar pair registered under POLYMER_TYPE: "DNA" must be recognised
 * via the DNA fallback in `get_na_natural_analog`.
 */
BOOST_AUTO_TEST_CASE(test_custom_dna_registered_sugar_pair_enables_toggle)
{
    constexpr std::string_view json = R"([
        {
            "symbol": "Yur",
            "polymer_type": "DNA",
            "natural_analog": "Yur",
            "smiles": "OC[C@H]1O[C@H]([H:1])[C@@H]([OH:2])[C@@H]1O",
            "core_smiles": "OC[C@H]1O[C@H]([H:1])[C@@H]([OH:2])[C@@H]1O",
            "name": "Custom Yur",
            "monomer_type": "Backbone",
            "author": "test",
            "pdbcode": "YUR"
        },
        {
            "symbol": "dYur",
            "polymer_type": "DNA",
            "natural_analog": "Yur",
            "smiles": "OC[C@H]1O[C@H]([H:1])C[C@@H]1[OH:2]",
            "core_smiles": "OC[C@H]1O[C@H]([H:1])C[C@@H]1[OH:2]",
            "name": "Custom dYur",
            "monomer_type": "Backbone",
            "author": "test",
            "pdbcode": "DYR"
        }
    ])";
    CustomDbScope scope(json);

    auto mol = rdkit_extensions::to_rdkit("RNA1{[Yur]}$$$$V2.0");
    BOOST_REQUIRE(get_monomer_type(mol->getAtomWithIdx(0)) ==
                  MonomerType::NA_SUGAR);

    TestMonomerContextMenu menu;
    MutateRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {},
                         mol->getAtomWithIdx(0));
    menu.updateActions();

    auto* toggle = find_action(menu, "Change Yur to dYur");
    BOOST_REQUIRE(toggle != nullptr);
    BOOST_TEST(toggle->isEnabled());
    toggle->trigger();

    BOOST_REQUIRE(capture.last_mutations.size() == 1u);
    BOOST_TEST(capture.last_mutations[0].helm_symbol == "dYur");
}

/**
 * "Add Complementary Sequence" hides for peptide selections — it only
 * makes sense for nucleic acids.
 */
BOOST_AUTO_TEST_CASE(test_complement_action_hidden_for_peptide_selection)
{
    auto mol = rdkit_extensions::to_rdkit("PEPTIDE1{A}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0}), {}, {}, {}, {});
    menu.updateActions();

    auto* action = find_action(menu, "Add Complementary Sequence");
    BOOST_REQUIRE(action != nullptr);
    BOOST_TEST(!action->isVisible());
}

/**
 * NA selection containing only sugars/phosphates (no bases) hides the
 * action — there's nothing to complement.
 */
BOOST_AUTO_TEST_CASE(test_complement_action_hidden_for_no_bases)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    // atoms 0=R sugar, 1=A base, 2=P phosphate
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {0, 2}), {}, {}, {}, {});
    menu.updateActions();

    auto* action = find_action(menu, "Add Complementary Sequence");
    BOOST_REQUIRE(action != nullptr);
    BOOST_TEST(!action->isVisible());
}

/**
 * A single NA base selection enables the action.
 */
BOOST_AUTO_TEST_CASE(test_complement_action_visible_for_single_base)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    menu.setContextItems(atoms_from(*mol, {1}), {}, {}, {}, {});
    menu.updateActions();

    auto* action = find_action(menu, "Add Complementary Sequence");
    BOOST_REQUIRE(action != nullptr);
    BOOST_TEST(action->isVisible());
    BOOST_TEST(action->isEnabled());
}

/**
 * Selection containing both bases and backbone monomers stays enabled —
 * the action filters down to the bases when it emits.
 */
BOOST_AUTO_TEST_CASE(test_complement_action_visible_for_mixed_na_selection)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}$$$$V2.0");
    TestMonomerContextMenu menu;
    ComplementRequestCapture capture;
    capture.connectTo(menu);
    menu.setContextItems(atoms_from(*mol, {0, 1, 2}), {}, {}, {}, {});
    menu.updateActions();

    auto* action = find_action(menu, "Add Complementary Sequence");
    BOOST_REQUIRE(action != nullptr);
    BOOST_TEST(action->isVisible());
    BOOST_TEST(action->isEnabled());
    action->trigger();

    BOOST_TEST(capture.call_count == 1);
    BOOST_REQUIRE(capture.last_bases.size() == 1u);
    BOOST_TEST(*capture.last_bases.begin() == mol->getAtomWithIdx(1));
}

/**
 * Mixing an NA base with a peptide hides the action — the spec restricts
 * it to "bases or full NAs only".
 */
BOOST_AUTO_TEST_CASE(test_complement_action_hidden_for_mixed_with_peptide)
{
    auto mol = rdkit_extensions::to_rdkit("RNA1{R(A)P}|PEPTIDE1{A}$$$$V2.0");
    TestMonomerContextMenu menu;
    // atom 1 = RNA base A, atom 3 = PEPTIDE A (positions depend on parse
    // order; verify with a quick monomer-type check).
    std::unordered_set<const RDKit::Atom*> sel;
    for (const auto* a : mol->atoms()) {
        if (get_monomer_type(a) == MonomerType::NA_BASE ||
            get_monomer_type(a) == MonomerType::PEPTIDE) {
            sel.insert(a);
        }
    }
    menu.setContextItems(sel, {}, {}, {}, {});
    menu.updateActions();

    auto* action = find_action(menu, "Add Complementary Sequence");
    BOOST_REQUIRE(action != nullptr);
    BOOST_TEST(!action->isVisible());
}

} // namespace sketcher
} // namespace schrodinger
