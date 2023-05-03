/* -------------------------------------------------------------------------
 * Tests class schrodinger::sketcher::CommonAtomPropertiesWidget
 --------------------------------------------------------------------------- */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE edit_atom_properties

#include <QPushButton>
#include <QSignalSpy>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "../test_sketcherScene.h"
#include "schrodinger/sketcher/Atom.h"
#include "schrodinger/sketcher/ChemicalKnowledge.h"
#include "schrodinger/sketcher/dialog/edit_atom_properties.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/ui/ui_common_atom_properties_widget.h"
#include "schrodinger/sketcher/ui/ui_edit_atom_properties.h"

using namespace schrodinger::sketcher;

BOOST_GLOBAL_FIXTURE(Test_Sketcher_global_fixture);

// Allow testing of equality operator
BOOST_TEST_DONT_PRINT_LOG_VALUE(EnhancedStereoType)
BOOST_TEST_DONT_PRINT_LOG_VALUE(QueryType)

/**
 * Subclass that allows access to protected data for testing purposes.
 */
class TestCommonAtomPropertiesWidget : public CommonAtomPropertiesWidget
{
  public:
    Ui::CommonAtomPropertiesWidget* getUI() const
    {
        return ui.get();
    }
    using CommonAtomPropertiesWidget::setEnhancedStereoTypeComboValue;
};

class TestEditAtomPropertiesDialog : public EditAtomPropertiesDialog
{
  public:
    TestEditAtomPropertiesDialog(SketcherModel* model, sketcherAtom& atom) :
        EditAtomPropertiesDialog(model, atom)
    {
    }
    Ui::CommonAtomPropertiesWidget* getAtomPropsUI() const
    {
        return ui->atom_common_props_wdg->ui.get();
    }
    Ui::EditAtomPropertiesDialog* getUI() const
    {
        return ui.get();
    }
    using EditAtomPropertiesDialog::accept;
    using EditAtomPropertiesDialog::getQueryTypeComboValue;
    using EditAtomPropertiesDialog::m_atom;
    using EditAtomPropertiesDialog::reset;
    using EditAtomPropertiesDialog::setQueryTypeComboValue;
    using EditAtomPropertiesDialog::updateOKButtonEnabled;
};

/**
 * Verify proper synchronization with model.
 */
BOOST_AUTO_TEST_CASE(test_CommonAtomPropertiesWidget)
{

    sketcherAtom atom;

    TestCommonAtomPropertiesWidget widget;
    auto ui = widget.getUI();

    // The clear button seems to leak in this test, so let's disable it.
    ui->isotope_le->setClearButtonEnabled(false);

    widget.readAtomInfo(atom);
    BOOST_TEST(ui->isotope_le->text().toStdString() == "");
    BOOST_TEST(ui->charge_sb->value() == 0);
    BOOST_TEST(ui->unpaired_elec_sb->value() == 0);
    BOOST_TEST(ui->enhanced_stereo_combo->currentText().toStdString() == "ABS");
    BOOST_TEST(ui->enhanced_stereo_sb->text().toStdString() == ""); // cleared
    BOOST_TEST(!ui->enhanced_stereo_sb->isEnabled());

    // Update UI elements
    ui->isotope_le->setText("13");
    ui->charge_sb->setValue(-2);
    ui->unpaired_elec_sb->setValue(1);
    ui->enhanced_stereo_combo->setCurrentIndex(2);
    BOOST_TEST(ui->enhanced_stereo_combo->currentText().toStdString() == "OR");
    ui->enhanced_stereo_sb->setValue(2);

    // Push data back out to the atom
    widget.writeAtomInfo(atom);
    BOOST_TEST(atom.getIsotope() == 13);
    BOOST_TEST(atom.getCharge() == -2);
    BOOST_TEST(atom.getUnpairedElectronsN() == 1);
    auto enh_st = atom.getEnhancedStereo();
    BOOST_TEST(enh_st.type == EnhancedStereoType::OR);
    BOOST_TEST(enh_st.group_id == 2);

    // Test setting isotope out of range
    ui->isotope_le->setText("1000");
    widget.writeAtomInfo(atom);
    BOOST_TEST(atom.getIsotope() == MAX_ISOTOPE_VALUE);
    ui->isotope_le->setText("-1000");
    widget.writeAtomInfo(atom);
    BOOST_TEST(atom.getIsotope() == MAX_ISOTOPE_VALUE);
}

/**
 * Verify that the widget properly reads/writes enhanced stereo information.
 */
BOOST_AUTO_TEST_CASE(read_and_write_atom_info)
{

    // We need to create an atom in the scene so that it will acquire the
    // necessary stereochemical properties
    testSketcherScene scene;
    scene.importText("C(F)(Cl)(Br)(I)");
    std::vector<sketcherAtom*> atoms;
    scene.getAtoms(atoms);
    sketcherAtom* c_atom = nullptr;
    for (auto atom : atoms) {
        if (atom->getAtomType() == static_cast<int>(Element::C)) {
            c_atom = atom;
            break;
        }
    }

    BOOST_TEST(c_atom->getNeighbors().size() == 4);
    c_atom->getBonds()[0]->setType(sketcherBond::BondType::SINGLE_UP);
    scene.forceUpdateStructure();
    BOOST_TEST(c_atom->getRDKAtom() != nullptr);

    // Verify that we can assign stereo properties on the atom itself
    EnhancedStereo exp_enhanced_stereo{EnhancedStereoType::AND, 1};
    c_atom->setEnhancedStereo(exp_enhanced_stereo);
    auto enhanced_stereo = c_atom->getEnhancedStereo();
    BOOST_TEST(enhanced_stereo.type == exp_enhanced_stereo.type);
    BOOST_TEST(enhanced_stereo.group_id == exp_enhanced_stereo.group_id);

    TestCommonAtomPropertiesWidget widget;
    widget.readAtomInfo(*c_atom);
    auto ui = widget.getUI();
    BOOST_TEST(ui->enhanced_stereo_combo->currentText().toStdString() == "AND");
    BOOST_TEST(ui->enhanced_stereo_sb->value() == exp_enhanced_stereo.group_id);

    widget.setEnhancedStereoTypeComboValue(EnhancedStereoType::ABS);
    widget.writeAtomInfo(*c_atom);

    enhanced_stereo = c_atom->getEnhancedStereo();
    BOOST_TEST(enhanced_stereo.type == EnhancedStereoType::ABS);
    BOOST_TEST(enhanced_stereo.group_id == 0);
}

/**
 * Verify that the widget properly resets inappropriate properties when
 * writing atom properties to the rdk atom.
 */
BOOST_AUTO_TEST_CASE(write_query_info)
{
    testSketcherScene scene;
    scene.importText("[CH+3]");
    auto atom = scene.quickGetAtoms().front();

    TestEditAtomPropertiesDialog dlg(scene.getModel(), *atom);
    auto ui = dlg.getUI();
    ui->set_as_query_rb->click();
    dlg.setQueryTypeComboValue(QueryType::RGROUP);
    dlg.accept();

    // That is, the charge is reset when turned into an rgroup atom
    BOOST_TEST(atom->isRGroup());
    BOOST_TEST(atom->getCharge() == 0);
}

/**
 * Verify that the dialog will not modify atoms if the user provided invalid
 * element or query information in a line edit.
 */
BOOST_AUTO_TEST_CASE(test_write_bad_info)
{
    testSketcherScene scene;
    scene.importText("C");
    auto atom = scene.quickGetAtoms().front();

    TestEditAtomPropertiesDialog dlg(scene.getModel(), *atom);
    auto ui = dlg.getUI();
    auto atom_props_ui = dlg.getAtomPropsUI();
    ui->set_as_query_rb->click();

    dlg.setQueryTypeComboValue(QueryType::SPECIFIC_ELEMENT);
    ui->specific_element_le->setText("whatever");
    atom_props_ui->charge_sb->setValue(2);
    dlg.accept();
    BOOST_TEST(atom->getAtomType() == static_cast<int>(Element::C));
    BOOST_TEST(atom->getCharge() == 0);

    dlg.setQueryTypeComboValue(QueryType::ALLOWED_LIST);
    ui->element_list_le->setText("whatever");
    BOOST_TEST(atom_props_ui->charge_sb->value() == 2);
    dlg.accept();
    BOOST_TEST(atom->getAtomType() == static_cast<int>(Element::C));
    BOOST_TEST(atom->getCharge() == 0);

    dlg.setQueryTypeComboValue(QueryType::NOT_ALLOWED_LIST);
    BOOST_TEST(ui->element_list_le->text() == QString("whatever"));
    BOOST_TEST(atom_props_ui->charge_sb->value() == 2);
    dlg.accept();
    BOOST_TEST(atom->getAtomType() == static_cast<int>(Element::C));
    BOOST_TEST(atom->getCharge() == 0);
}

void check_atom_query_page(TestEditAtomPropertiesDialog& dlg)
{
    auto query_type = dlg.getQueryTypeComboValue();
    auto ui = dlg.getUI();
    ui->atom_query_stacked_wdg->setCurrentWidget(ui->edit_query_page);

    auto in_specific_element = query_type == QueryType::SPECIFIC_ELEMENT;
    auto in_element_list = query_type == QueryType::ALLOWED_LIST ||
                           query_type == QueryType::NOT_ALLOWED_LIST;
    auto in_wildcard = query_type == QueryType::WILDCARD;
    auto in_rgroup = query_type == QueryType::RGROUP;

    auto symbol =
        QString(atomic_number_to_symbol(dlg.m_atom.getAtomicNumber()).c_str());

    BOOST_TEST(ui->periodic_table_btn->isVisibleTo(&dlg) ==
               (in_specific_element || in_element_list));
    BOOST_TEST(ui->wildcard_combo->isVisibleTo(&dlg) == in_wildcard);
    BOOST_TEST(ui->specific_element_le->isVisibleTo(&dlg) ==
               in_specific_element);
    BOOST_TEST(ui->element_list_le->isVisibleTo(&dlg) == in_element_list);
    BOOST_TEST(ui->rgroup_sb->isVisibleTo(&dlg) == in_rgroup);
    BOOST_TEST(ui->specific_element_le->placeholderText() ==
               QString("(element)"));
    BOOST_TEST(ui->element_list_le->placeholderText() ==
               QString("(comma-separated list)"));
    BOOST_TEST(ui->specific_element_le->text() == symbol);
    BOOST_TEST(ui->element_list_le->text() == symbol);
}

/**
 * Verify appearance and behavior for the atom query page.
 */
BOOST_AUTO_TEST_CASE(atom_query_page)
{

    sketcherAtom atom;
    atom.setAtomType(6);
    SketcherModel model;
    TestEditAtomPropertiesDialog dlg(&model, atom);
    BOOST_TEST(dlg.getQueryTypeComboValue() ==
               QueryType::SPECIFIC_ELEMENT); // Default
    for (auto query_type :
         {QueryType::WILDCARD, QueryType::ALLOWED_LIST,
          QueryType::SPECIFIC_ELEMENT, QueryType::NOT_ALLOWED_LIST}) {
        dlg.setQueryTypeComboValue(query_type);
        BOOST_TEST(dlg.getQueryTypeComboValue() == query_type);
        check_atom_query_page(dlg);
    }
}

/** Test if dialog fires structure changed signal on pressing OK
 * Method is intended to verify only signal.
 */
BOOST_AUTO_TEST_CASE(edit_atom_property_structure_change_signal)
{

    testSketcherScene scene;
    scene.importText("CCCCCC");
    std::vector<sketcherAtom*> atoms;
    scene.getAtoms(atoms);
    auto atom = atoms[0];
    TestEditAtomPropertiesDialog dlg(scene.getModel(), *atom);
    auto ui = dlg.getUI();
    ui->element_le->setText("Br");
    QSignalSpy structure_change_signal_spy(&scene,
                                           &sketcherScene::structureChanged);
    dlg.accept();
    BOOST_CHECK_EQUAL(structure_change_signal_spy.count(), 1);
}

/**
 * Verify behavior of R group spinbox
 */
BOOST_AUTO_TEST_CASE(rgroup_spinbox)
{

    testSketcherScene scene;
    scene.importText("CC");
    std::vector<sketcherAtom*> atoms;
    scene.getAtoms(atoms);
    auto element_atom = atoms[0];
    auto rgroup_atom = atoms[1];
    rgroup_atom->setAsRGroup(1);

    // When a non-R group is imported, the spin box should show the next
    // available R group number from the scene
    TestEditAtomPropertiesDialog element_dlg(scene.getModel(), *element_atom);
    auto ui = element_dlg.getUI();
    BOOST_TEST(ui->rgroup_sb->value() == 2);
    ui->atom_query_stacked_wdg->setCurrentWidget(ui->edit_query_page);
    element_dlg.setQueryTypeComboValue(QueryType::RGROUP);
    element_dlg.accept();
    BOOST_TEST(element_atom->isRGroup());
    BOOST_TEST(element_atom->getRGroupNumber() == 2);

    // When an R group is imported, the spin box should show the R group number
    // that is currently used by the atom object
    TestEditAtomPropertiesDialog rgroup_dlg(scene.getModel(), *rgroup_atom);
    BOOST_TEST(rgroup_dlg.getUI()->rgroup_sb->value() == 1);
}

BOOST_AUTO_TEST_CASE(test_ring_count_spinboxes)
{
    testSketcherScene scene;
    scene.importText("CC");
    auto atom = scene.quickGetAtoms().front();

    // When an R group is imported, the spin box should show the R group number
    // that is currently used by the atom object
    TestEditAtomPropertiesDialog dlg(scene.getModel(), *atom);
    auto ui = dlg.getUI();

    auto set_text = [](QComboBox* combo, const QString& text) {
        combo->setCurrentIndex(combo->findText(text));
    };

    for (const auto& pair : std::vector<std::pair<QComboBox*, QSpinBox*>>{
             {ui->ring_count_combo, ui->ring_count_sb},
             {ui->ring_bond_count_combo, ui->ring_bond_count_sb}}) {
        auto combo = pair.first;
        auto sb = pair.second;
        // combo initially set to (any)
        BOOST_TEST(!sb->isEnabled());
        set_text(combo, "exactly");
        BOOST_TEST(sb->isEnabled());
        BOOST_TEST(sb->value() == 0);
        sb->setValue(4);
        BOOST_TEST(sb->value() == 4);
        set_text(combo, "(any)");
        BOOST_TEST(!sb->isEnabled());
        set_text(combo, "exactly");
        BOOST_TEST(sb->isEnabled());
        BOOST_TEST(sb->value() == 4);
    }
}

BOOST_AUTO_TEST_CASE(test_updateOKButtonEnabled)
{
    testSketcherScene scene;
    scene.importText("C");
    auto atom = scene.quickGetAtoms().front();
    atom->setAsRGroup(1);
    TestEditAtomPropertiesDialog dlg(scene.getModel(), *atom);
    auto ui = dlg.getUI();
    auto ok_button = ui->buttonBox->button(QDialogButtonBox::Ok);
    // Query mode tests
    ui->set_as_query_rb->click();
    // Test RGROUP query type.
    dlg.setQueryTypeComboValue(QueryType::RGROUP);
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), true);

    // Test WILDCARD query type
    dlg.setQueryTypeComboValue(QueryType::WILDCARD);
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), true);

    // Test ALLOWED_LIST query type.
    dlg.setQueryTypeComboValue(QueryType::ALLOWED_LIST);
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), false);
    ui->element_list_le->setText("Br");
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), true);

    // Test NOT_ALLOWED_LIST query type.
    dlg.setQueryTypeComboValue(QueryType::NOT_ALLOWED_LIST);
    BOOST_CHECK_EQUAL(ok_button->isEnabled(),
                      !ui->element_list_le->text().isEmpty());
    ui->element_list_le->setText("Br");
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), true);
    ui->element_list_le->setText("");
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), false);

    // Test SPECIFIC_ELEMENT query type.
    dlg.setQueryTypeComboValue(QueryType::SPECIFIC_ELEMENT);
    BOOST_CHECK_EQUAL(ok_button->isEnabled(),
                      !ui->specific_element_le->text().isEmpty());
    ui->specific_element_le->setText("Br");
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), true);
    ui->specific_element_le->setText("");
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), false);

    // Atom mode tests
    ui->set_as_atom_rb->click();
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), true);
    ui->element_le->setText("");
    BOOST_CHECK_EQUAL(ok_button->isEnabled(), false);
}

BOOST_AUTO_TEST_CASE(test_specific_element_advanced_query)
{
    testSketcherScene scene;
    scene.importText("C");
    auto atom = scene.quickGetAtoms().front();

    TestEditAtomPropertiesDialog dlg(scene.getModel(), *atom);
    auto ui = dlg.getUI();
    BOOST_TEST(ui->atom_query_stacked_wdg->currentWidget() ==
               ui->edit_atom_page);

    ui->set_as_query_rb->click();
    BOOST_TEST(dlg.getQueryTypeComboValue() == QueryType::SPECIFIC_ELEMENT);
    BOOST_TEST(ui->specific_element_le->text().toStdString() == "C");
    ui->total_H_le->setText("1");
    dlg.accept();

    dlg.reset();
    BOOST_TEST(ui->atom_query_stacked_wdg->currentWidget() ==
               ui->edit_query_page);
    BOOST_TEST(dlg.getQueryTypeComboValue() == QueryType::SPECIFIC_ELEMENT);
    BOOST_TEST(ui->specific_element_le->text().toStdString() == "C");
}

BOOST_AUTO_TEST_CASE(test_elementList)
{
    // Tests the element list line edit by default is set to the selected atom
    testSketcherScene scene;
    scene.importText("C");
    auto atom = scene.quickGetAtoms().front();

    TestEditAtomPropertiesDialog dlg(scene.getModel(), *atom);
    auto ui = dlg.getUI();
    BOOST_TEST(ui->element_list_le->text().toStdString() == "C");
}