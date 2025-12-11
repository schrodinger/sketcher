
#define BOOST_TEST_MODULE test_edit_atom_properties

#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/dialog/edit_atom_properties.h"
#include "schrodinger/sketcher/ui/ui_edit_atom_properties.h"
#include "schrodinger/sketcher/widget/widget_utils.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);
BOOST_TEST_DONT_PRINT_LOG_VALUE(std::optional<int>)
BOOST_TEST_DONT_PRINT_LOG_VALUE(std::optional<unsigned int>)
BOOST_TEST_DONT_PRINT_LOG_VALUE(std::nullopt_t)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::QueryType)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::QueryAromaticity)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::QueryCount)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::EnhancedStereo)
BOOST_TEST_DONT_PRINT_LOG_VALUE(
    std::optional<schrodinger::sketcher::EnhancedStereo>)

namespace schrodinger
{
namespace sketcher
{

class TestDialog : public EditAtomPropertiesDialog
{
  public:
    using EditAtomPropertiesDialog::EditAtomPropertiesDialog;
    using EditAtomPropertiesDialog::getDialogSettings;
    using EditAtomPropertiesDialog::ui;
};

/**
 * Make sure that loading a non-query atom into the dialog updates the widgets
 * to the correct property values
 */
BOOST_AUTO_TEST_CASE(test_load_atom)
{
    auto undo_stack = QUndoStack();
    auto mol_model = MolModel(&undo_stack);
    import_mol_text(&mol_model, "[N++]");
    auto* atom = mol_model.getMol()->getAtomWithIdx(0);
    auto dialog = TestDialog(atom, &mol_model);
    BOOST_TEST(dialog.ui->set_as_atom_rb->isChecked());
    BOOST_TEST(dialog.ui->atom_query_stacked_wdg->currentWidget() ==
               dialog.ui->edit_atom_page);
    BOOST_TEST(dialog.ui->atom_isotope_sb->optionalValue() == std::nullopt);
    BOOST_TEST(dialog.ui->atom_element_le->text() == "N");
    BOOST_TEST(dialog.ui->atom_charge_sb->value() == 2);
    BOOST_TEST(dialog.ui->atom_unpaired_sb->value() == 3);
}

/**
 * Make sure that loading a query atom into the dialog updates the widgets
 * to the correct property values
 */
BOOST_AUTO_TEST_CASE(test_load_query)
{
    auto undo_stack = QUndoStack();
    auto mol_model = MolModel(&undo_stack);
    import_mol_text(&mol_model, "[nR3++]");
    auto* atom = mol_model.getMol()->getAtomWithIdx(0);
    auto dialog = TestDialog(atom, &mol_model);
    BOOST_TEST(dialog.ui->set_as_query_rb->isChecked());
    BOOST_TEST(dialog.ui->atom_query_stacked_wdg->currentWidget() ==
               dialog.ui->edit_query_page);
    BOOST_TEST(dialog.ui->edit_query_tab_wdg->currentWidget() ==
               dialog.ui->general_tab);
    BOOST_TEST(dialog.ui->query_type_combo->currentData().value<QueryType>() ==
               QueryType::SPECIFIC_ELEMENT);
    BOOST_TEST(dialog.ui->query_element_le->text() == "N");
    BOOST_TEST(dialog.ui->query_isotope_sb->optionalValue() == std::nullopt);
    BOOST_TEST(dialog.ui->query_charge_sb->optionalValue() == 2);
    BOOST_TEST(dialog.ui->query_unpaired_sb->optionalValue() == std::nullopt);
    BOOST_TEST(dialog.ui->query_unpaired_sb->optionalValue() == std::nullopt);
    BOOST_TEST(
        dialog.ui->aromaticity_combo->currentData().value<QueryAromaticity>() ==
        QueryAromaticity::AROMATIC);
    BOOST_TEST(dialog.ui->ring_count_combo->currentData().value<QueryCount>() ==
               QueryCount::EXACTLY);
    BOOST_TEST(dialog.ui->ring_count_sb->value() == 3);
    BOOST_TEST(
        dialog.ui->ring_bond_count_combo->currentData().value<QueryCount>() ==
        QueryCount::ANY);
}

/**
 * Make sure that getDialogSettings returns an atom properties object with the
 * correct values
 */
BOOST_AUTO_TEST_CASE(test_getDialogSettings)
{
    auto undo_stack = QUndoStack();
    auto mol_model = MolModel(&undo_stack);
    import_mol_text(&mol_model, "C");
    auto* atom = mol_model.getMol()->getAtomWithIdx(0);
    auto dialog = TestDialog(atom, &mol_model);
    auto props = dialog.getDialogSettings();
    BOOST_TEST(!props->isQuery());
    auto* atom_props = static_cast<AtomProperties*>(props.get());
    BOOST_TEST(atom_props->element == Element::C);
    BOOST_TEST(atom_props->isotope == std::nullopt);
    BOOST_TEST(atom_props->charge == 0);

    dialog.ui->atom_element_le->setText("Mg");
    dialog.ui->atom_charge_sb->setValue(2);
    props = dialog.getDialogSettings();
    BOOST_TEST(!props->isQuery());
    atom_props = static_cast<AtomProperties*>(props.get());
    BOOST_TEST(atom_props->element == Element::MG);
    BOOST_TEST(atom_props->isotope == std::nullopt);
    BOOST_TEST(atom_props->charge == 2);

    dialog.ui->set_as_query_rb->click();
    set_combo_box_data(dialog.ui->query_type_combo,
                       QueryType::SPECIFIC_ELEMENT);
    dialog.ui->query_element_le->setText("N");
    dialog.ui->query_charge_sb->setOptionalValue(1);
    dialog.ui->query_unpaired_sb->setOptionalValue(2);
    set_combo_box_data(dialog.ui->ring_count_combo, QueryCount::EXACTLY);
    dialog.ui->ring_count_sb->setValue(3);
    props = dialog.getDialogSettings();
    BOOST_TEST(props->isQuery());
    auto* query_props = static_cast<AtomQueryProperties*>(props.get());
    BOOST_TEST(query_props->query_type == QueryType::SPECIFIC_ELEMENT);
    BOOST_TEST(query_props->element == Element::N);
    BOOST_TEST(query_props->charge == 1);
    BOOST_TEST(query_props->unpaired_electrons == 2);
    BOOST_TEST(query_props->ring_count_type == QueryCount::EXACTLY);
    BOOST_TEST(query_props->ring_count_exact_val == 3);

    // make sure that getDialogSettings ignores disabled widgets even if they
    // have a value in them - switching to "Allowed List" will disable the
    // "Unpaired electrons" spin box
    set_combo_box_data(dialog.ui->query_type_combo, QueryType::ALLOWED_LIST);
    props = dialog.getDialogSettings();
    BOOST_TEST(props->isQuery());
    query_props = static_cast<AtomQueryProperties*>(props.get());
    BOOST_TEST(query_props->query_type == QueryType::ALLOWED_LIST);
    BOOST_TEST(query_props->allowed_list == std::vector<Element>{Element::N});
    BOOST_TEST(query_props->charge == 1);
    BOOST_TEST(query_props->unpaired_electrons == std::nullopt);
    BOOST_TEST(query_props->ring_count_type == QueryCount::EXACTLY);
    BOOST_TEST(query_props->ring_count_exact_val == 3);

    // make sure that getDialog settings ignores hidden widgets, such as the
    // ring (bond) count spin box when the ring (bond) count combo box is set to
    // any or positive
    set_combo_box_data(dialog.ui->ring_bond_count_combo, QueryCount::EXACTLY);
    dialog.ui->ring_bond_count_sb->setValue(2);
    set_combo_box_data(dialog.ui->ring_bond_count_combo, QueryCount::ANY);
    set_combo_box_data(dialog.ui->ring_count_combo, QueryCount::POSITIVE);
    props = dialog.getDialogSettings();
    BOOST_TEST(props->isQuery());
    query_props = static_cast<AtomQueryProperties*>(props.get());
    BOOST_TEST(query_props->query_type == QueryType::ALLOWED_LIST);
    BOOST_TEST(query_props->allowed_list == std::vector<Element>{Element::N});
    BOOST_TEST(query_props->charge == 1);
    BOOST_TEST(query_props->unpaired_electrons == std::nullopt);
    BOOST_TEST(query_props->ring_count_type == QueryCount::POSITIVE);
    BOOST_TEST(query_props->ring_count_exact_val == 0);
    BOOST_TEST(query_props->ring_bond_count_type == QueryCount::ANY);
    BOOST_TEST(query_props->ring_bond_count_exact_val == 0);
}

/**
 * Make sure that getDialogSettings returns an atom properties object with the
 * correct enhanced stereo values
 */
BOOST_AUTO_TEST_CASE(test_getDialogSettings_enhanced_stereo)
{
    auto undo_stack = QUndoStack();
    auto mol_model = MolModel(&undo_stack);
    import_mol_text(&mol_model, "N[C@H](C)C(=O)O |&1:1|");
    auto* atom = mol_model.getMol()->getAtomWithIdx(1);
    auto dialog = TestDialog(atom, &mol_model);
    auto props = dialog.getDialogSettings();
    BOOST_TEST(!props->isQuery());
    auto* atom_props = static_cast<AtomProperties*>(props.get());
    BOOST_TEST(atom_props->element == Element::C);
    BOOST_TEST(atom_props->enhanced_stereo ==
               EnhancedStereo(RDKit::StereoGroupType::STEREO_AND, 1));

    set_combo_box_data(dialog.ui->atom_stereo_combo,
                       RDKit::StereoGroupType::STEREO_OR);
    dialog.ui->atom_stereo_sb->setValue(3);
    props = dialog.getDialogSettings();
    BOOST_TEST(!props->isQuery());
    atom_props = static_cast<AtomProperties*>(props.get());
    BOOST_TEST(atom_props->element == Element::C);
    BOOST_TEST(atom_props->enhanced_stereo ==
               EnhancedStereo(RDKit::StereoGroupType::STEREO_OR, 3));

    set_combo_box_data(dialog.ui->atom_stereo_combo,
                       RDKit::StereoGroupType::STEREO_ABSOLUTE);
    props = dialog.getDialogSettings();
    BOOST_TEST(!props->isQuery());
    atom_props = static_cast<AtomProperties*>(props.get());
    BOOST_TEST(atom_props->element == Element::C);
    BOOST_TEST(atom_props->enhanced_stereo ==
               EnhancedStereo(RDKit::StereoGroupType::STEREO_ABSOLUTE, 0));
}

/**
 * Make sure that settings are transferred between the atom page and the query
 * page when the page is switched. Also make sure that elements are transferred
 * between the query allowed list and the query specific element when we change
 * query type.
 */
BOOST_AUTO_TEST_CASE(test_transferring_settings)
{
    auto undo_stack = QUndoStack();
    auto mol_model = MolModel(&undo_stack);
    import_mol_text(&mol_model, "C");
    auto* atom = mol_model.getMol()->getAtomWithIdx(0);
    auto dialog = TestDialog(atom, &mol_model);
    dialog.ui->atom_charge_sb->setValue(2);

    dialog.ui->set_as_query_rb->click();
    BOOST_TEST(dialog.ui->query_type_combo->currentData().value<QueryType>() ==
               QueryType::SPECIFIC_ELEMENT);
    BOOST_TEST(dialog.ui->query_element_le->text() == "C");
    BOOST_TEST(dialog.ui->query_charge_sb->optionalValue() == 2);
    dialog.ui->query_charge_sb->setValue(1);

    dialog.ui->set_as_atom_rb->click();
    BOOST_TEST(dialog.ui->atom_charge_sb->value() == 1);
    BOOST_TEST(dialog.ui->atom_element_le->text() == "C");
    dialog.ui->atom_element_le->setText("N");

    dialog.ui->set_as_query_rb->click();
    set_combo_box_data(dialog.ui->query_type_combo, QueryType::ALLOWED_LIST);
    BOOST_TEST(dialog.ui->query_type_combo->currentData().value<QueryType>() ==
               QueryType::ALLOWED_LIST);
    BOOST_TEST(dialog.ui->element_list_le->text() == "N");
    dialog.ui->element_list_le->setText("Mg, P, K");

    // only the first element of the list should get transferred to the atom
    // element
    dialog.ui->set_as_atom_rb->click();
    BOOST_TEST(dialog.ui->atom_element_le->text() == "Mg");

    // we shouldn't truncate the query element list since the atom element is
    // the same as the first element of the allowed list
    dialog.ui->set_as_query_rb->click();
    BOOST_TEST(dialog.ui->element_list_le->text() == "Mg, P, K");

    dialog.ui->set_as_atom_rb->click();
    BOOST_TEST(dialog.ui->atom_element_le->text() == "Mg");
    dialog.ui->atom_element_le->setText("O");

    // this time, the element list should be completely replaced by the atom
    // element since it was changed
    dialog.ui->set_as_query_rb->click();
    BOOST_TEST(dialog.ui->element_list_le->text() == "O");

    // check that transferring settings between the query element list and the
    // query specific element works
    set_combo_box_data(dialog.ui->query_type_combo,
                       QueryType::SPECIFIC_ELEMENT);
    BOOST_TEST(dialog.ui->query_element_le->text() == "O");
    dialog.ui->query_element_le->setText("C");
    set_combo_box_data(dialog.ui->query_type_combo, QueryType::ALLOWED_LIST);
    BOOST_TEST(dialog.ui->element_list_le->text() == "C");
    dialog.ui->element_list_le->setText("F, Cl, Br");
    set_combo_box_data(dialog.ui->query_type_combo,
                       QueryType::SPECIFIC_ELEMENT);
    BOOST_TEST(dialog.ui->query_element_le->text() == "F");
    set_combo_box_data(dialog.ui->query_type_combo, QueryType::ALLOWED_LIST);
    BOOST_TEST(dialog.ui->element_list_le->text() == "F, Cl, Br");
    set_combo_box_data(dialog.ui->query_type_combo,
                       QueryType::SPECIFIC_ELEMENT);
    BOOST_TEST(dialog.ui->query_element_le->text() == "F");
    dialog.ui->query_element_le->setText("Ne");
    set_combo_box_data(dialog.ui->query_type_combo, QueryType::ALLOWED_LIST);
    BOOST_TEST(dialog.ui->element_list_le->text() == "Ne");
}

} // namespace sketcher
} // namespace schrodinger
