#define BOOST_TEST_MODULE Test_Sketcher

#include <QDialogButtonBox>
#include <QPushButton>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/custom_monomer_dialog.h"
#include "schrodinger/sketcher/sketcher_widget.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

static QPushButton* get_ok_button(const QWidget& dialog)
{
    auto* button_box = dialog.findChild<QDialogButtonBox*>("button_box");
    BOOST_REQUIRE(button_box != nullptr);
    return button_box->button(QDialogButtonBox::Ok);
}

/**
 * Make sure that the OK button is enabled only when the dialog is non-empty
 */
BOOST_AUTO_TEST_CASE(custom_monomer_dialog_validation_and_acceptance)
{
    CustomMonomerDialog dialog(rdkit_extensions::ChainType::PEPTIDE);
    auto* sketcher = dialog.findChild<SketcherWidget*>();
    auto* ok_button = get_ok_button(dialog);
    BOOST_REQUIRE(sketcher != nullptr);

    BOOST_TEST(dialog.windowTitle() == "Sketch Custom Peptide Monomer");
    BOOST_TEST(!ok_button->isEnabled());

    dialog.addSMILES("CC");
    BOOST_TEST(ok_button->isEnabled());

    bool monomer_accepted = false;
    auto accepted_chain_type = rdkit_extensions::ChainType::CHEM;
    QObject::connect(
        &dialog, &CustomMonomerDialog::customMonomerAccepted, &dialog,
        [&monomer_accepted, &accepted_chain_type](
            const std::string&, const rdkit_extensions::ChainType type) {
            monomer_accepted = true;
            accepted_chain_type = type;
        });
    dialog.accept();
    BOOST_TEST(monomer_accepted);
    BOOST_TEST(static_cast<int>(accepted_chain_type) ==
               static_cast<int>(rdkit_extensions::ChainType::PEPTIDE));
}

BOOST_AUTO_TEST_CASE(custom_monomer_dialog_titles_reflect_chain_type)
{
    CustomMonomerDialog peptide_dialog(rdkit_extensions::ChainType::PEPTIDE);
    CustomMonomerDialog nucleic_acid_dialog(rdkit_extensions::ChainType::RNA);
    CustomMonomerDialog chem_dialog(rdkit_extensions::ChainType::CHEM);

    BOOST_TEST(peptide_dialog.windowTitle() == "Sketch Custom Peptide Monomer");
    BOOST_TEST(nucleic_acid_dialog.windowTitle() ==
               "Sketch Custom Nucleic Acid Monomer");
    BOOST_TEST(chem_dialog.windowTitle() == "Sketch Custom Chem Monomer");
}

} // namespace sketcher
} // namespace schrodinger
