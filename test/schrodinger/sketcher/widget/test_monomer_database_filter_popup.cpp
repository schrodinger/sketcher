#define BOOST_TEST_MODULE test_monomer_database_filter_popup
#include <boost/test/unit_test.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include <QAbstractItemView>
#include <QApplication>
#include <QComboBox>
#include <QFont>
#include <QLabel>
#include <QTest>
#include <QToolButton>

#include "../test_common.h"
#include "schrodinger/sketcher/model/monomer_database_filter.h"
#include "schrodinger/sketcher/widget/monomer_database_filter_popup.h"
#include "schrodinger/sketcher/widget/multi_select_combo_box.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

static const auto PEPTIDE = MonomerDbPolymerType::PEPTIDE;
static const auto NUCLEIC_ACID = MonomerDbPolymerType::NUCLEIC_ACID;
static const auto CHEM = MonomerDbPolymerType::CHEM;

/**
 * A filter popup that's open next to an anchor widget
 */
struct FilterPopupFixture {
    FilterPopupFixture(
        const MonomerDatabaseFilter& filter = MonomerDatabaseFilter()) :
        popup(&anchor)
    {
        anchor.resize(30, 30);
        anchor.show();
        popup.showPopup(filter);
        BOOST_REQUIRE(popup.isVisible());
        polymer_type_combo =
            popup.findChild<MultiSelectComboBox*>("polymer_type_combo");
        monomer_type_combo = popup.findChild<QComboBox*>("monomer_type_combo");
        natural_analog_combo =
            popup.findChild<QComboBox*>("natural_analog_combo");
        filters_applied_lbl = popup.findChild<QLabel*>("filters_applied_lbl");
        BOOST_REQUIRE(polymer_type_combo != nullptr);
        BOOST_REQUIRE(monomer_type_combo != nullptr);
        BOOST_REQUIRE(natural_analog_combo != nullptr);
        BOOST_REQUIRE(filters_applied_lbl != nullptr);
    }

    ~FilterPopupFixture()
    {
        popup.hide();
    }

    /**
     * Check or uncheck a polymer type as the user would
     */
    void setPolymerTypeChecked(const QString& text, bool checked)
    {
        int row = polymer_type_combo->findText(text);
        BOOST_REQUIRE(row >= 0);
        polymer_type_combo->setItemChecked(row, checked);
    }

    /**
     * Choose an item in a single-select combo box as the user would
     */
    void choose(QComboBox* combo, const QString& text)
    {
        int row = combo->findText(text);
        BOOST_REQUIRE_MESSAGE(row >= 0, text.toStdString());
        combo->setCurrentIndex(row);
    }

    std::vector<std::string> getItemTexts(QComboBox* combo)
    {
        std::vector<std::string> texts;
        for (int row = 0; row < combo->count(); ++row) {
            texts.push_back(combo->itemText(row).toStdString());
        }
        return texts;
    }

    /**
     * Check that the combo box's dropdown draws none of its items in bold,
     * regardless of whether the combo box itself is bold. The dropdown draws
     * an item with the item's Qt::FontRole data if it has any, and with the
     * combo box's font otherwise.
     */
    void checkItemsNotBold(QComboBox* combo)
    {
        for (int row = 0; row < combo->count(); ++row) {
            const auto font_data = combo->itemData(row, Qt::FontRole);
            BOOST_TEST_REQUIRE(font_data.isValid());
            BOOST_TEST(!font_data.value<QFont>().bold());
        }
    }

    /**
     * Check what each combo box shows, whether it's enabled and whether it's
     * bold. The combo boxes' dropdown items are never bold.
     */
    void checkCombos(const QString& polymer_types, bool polymer_types_bold,
                     const QString& monomer_type, bool monomer_type_enabled,
                     bool monomer_type_bold, const QString& analog,
                     bool analog_enabled, bool analog_bold)
    {
        BOOST_TEST(polymer_type_combo->getDisplayText() == polymer_types);
        BOOST_TEST(polymer_type_combo->font().bold() == polymer_types_bold);
        BOOST_TEST(monomer_type_combo->currentText() == monomer_type);
        BOOST_TEST(monomer_type_combo->isEnabled() == monomer_type_enabled);
        BOOST_TEST(monomer_type_combo->font().bold() == monomer_type_bold);
        BOOST_TEST(natural_analog_combo->currentText() == analog);
        BOOST_TEST(natural_analog_combo->isEnabled() == analog_enabled);
        BOOST_TEST(natural_analog_combo->font().bold() == analog_bold);
        for (auto* combo : {static_cast<QComboBox*>(polymer_type_combo),
                            monomer_type_combo, natural_analog_combo}) {
            checkItemsNotBold(combo);
        }
    }

    QWidget anchor;
    MonomerDatabaseFilterPopup popup;
    MultiSelectComboBox* polymer_type_combo = nullptr;
    QComboBox* monomer_type_combo = nullptr;
    QComboBox* natural_analog_combo = nullptr;
    QLabel* filters_applied_lbl = nullptr;
};

/**
 * The combo boxes list the expected items, and the natural analogs depend on
 * the polymer type that the natural analog applies to
 */
BOOST_AUTO_TEST_CASE(combo_contents)
{
    FilterPopupFixture fixture{MonomerDatabaseFilter(PEPTIDE | CHEM)};
    BOOST_TEST(
        fixture.getItemTexts(fixture.polymer_type_combo) ==
            (std::vector<std::string>{"Peptide", "Nucleic acid", "Structure"}),
        boost::test_tools::per_element());
    BOOST_TEST(fixture.getItemTexts(fixture.monomer_type_combo) ==
                   (std::vector<std::string>{"Any", "Backbone", "Terminal"}),
               boost::test_tools::per_element());
    auto analogs = fixture.getItemTexts(fixture.natural_analog_combo);
    BOOST_REQUIRE(analogs.size() == 23);
    BOOST_TEST(
        (std::vector<std::string>(analogs.begin(), analogs.begin() + 3)) ==
            (std::vector<std::string>{"Any or None", "Any", "None"}),
        boost::test_tools::per_element());
    // the 20 standard amino acids, including N (asparagine) but not X
    BOOST_TEST(std::count(analogs.begin(), analogs.end(), "N") == 1);
    BOOST_TEST(std::count(analogs.begin(), analogs.end(), "X") == 0);
    BOOST_TEST(std::count(analogs.begin(), analogs.end(), "C") == 1);

    fixture.setPolymerTypeChecked("Nucleic acid", true);
    BOOST_TEST(fixture.getItemTexts(fixture.natural_analog_combo) ==
                   (std::vector<std::string>{"Any or None", "Any", "None"}),
               boost::test_tools::per_element());
    fixture.setPolymerTypeChecked("Peptide", false);
    BOOST_TEST(fixture.getItemTexts(fixture.natural_analog_combo) ==
                   (std::vector<std::string>{"Any or None", "Any", "None", "R",
                                             "P", "A", "C", "G", "U", "T"}),
               boost::test_tools::per_element());
}

/**
 * Monomer type is enabled only for {Peptide} and {Peptide, Structure}, and
 * Natural analog whenever exactly one of Peptide and Nucleic acid is checked
 */
BOOST_AUTO_TEST_CASE(enabling)
{
    struct Case {
        MonomerDbPolymerType polymer_types;
        bool monomer_type_enabled;
        bool analog_enabled;
    };
    for (const auto& [polymer_types, monomer_type_enabled, analog_enabled] :
         {Case{PEPTIDE, true, true}, Case{PEPTIDE | CHEM, true, true},
          Case{NUCLEIC_ACID, false, true},
          Case{NUCLEIC_ACID | CHEM, false, true}, Case{CHEM, false, false},
          Case{PEPTIDE | NUCLEIC_ACID, false, false},
          Case{MonomerDbPolymerType::ALL, false, false}}) {
        BOOST_TEST_CONTEXT(
            "polymer types=" << static_cast<unsigned>(polymer_types))
        {
            FilterPopupFixture fixture{MonomerDatabaseFilter(polymer_types)};
            BOOST_TEST(fixture.monomer_type_combo->isEnabled() ==
                       monomer_type_enabled);
            BOOST_TEST(fixture.natural_analog_combo->isEnabled() ==
                       analog_enabled);
        }
    }
}

/**
 * The page defaults count one filter with Polymer type in bold, and "Any" is
 * displayed when all three types are checked
 */
BOOST_AUTO_TEST_CASE(bold_and_count)
{
    FilterPopupFixture fixture{MonomerDatabaseFilter(PEPTIDE | CHEM)};
    fixture.checkCombos("Peptide, Structure", true, "Any", true, false,
                        "Any or None", true, false);
    BOOST_TEST(fixture.filters_applied_lbl->text() == "1 filter applied");

    fixture.choose(fixture.monomer_type_combo, "Backbone");
    fixture.choose(fixture.natural_analog_combo, "C");
    fixture.checkCombos("Peptide, Structure", true, "Backbone", true, true, "C",
                        true, true);
    BOOST_TEST(fixture.filters_applied_lbl->text() == "3 filters applied");

    // checking every polymer type makes both other combos inapplicable
    fixture.setPolymerTypeChecked("Nucleic acid", true);
    fixture.checkCombos("Any", false, "Any", false, false, "Any or None", false,
                        false);
    BOOST_TEST(fixture.filters_applied_lbl->text() == "No filters applied");
}

/**
 * Every user change emits the filter, and Clear All produces the most
 * permissive filter
 */
BOOST_AUTO_TEST_CASE(filter_changed_and_clear_all)
{
    FilterPopupFixture fixture{MonomerDatabaseFilter(PEPTIDE)};
    std::vector<MonomerDatabaseFilter> emitted;
    QObject::connect(&fixture.popup, &MonomerDatabaseFilterPopup::filterChanged,
                     [&emitted](const MonomerDatabaseFilter& filter) {
                         emitted.push_back(filter);
                     });
    fixture.choose(fixture.monomer_type_combo, "Terminal");
    BOOST_REQUIRE(emitted.size() == 1);
    BOOST_CHECK(emitted.back().getMonomerType() == MonomerTypeFilter::TERMINAL);
    fixture.choose(fixture.natural_analog_combo, "None");
    BOOST_REQUIRE(emitted.size() == 2);
    BOOST_CHECK(emitted.back().getNaturalAnalog() == NaturalAnalogFilter::NONE);
    BOOST_CHECK(fixture.popup.getFilter() == emitted.back());
    fixture.setPolymerTypeChecked("Structure", true);
    BOOST_REQUIRE(emitted.size() == 3);
    BOOST_CHECK(emitted.back().getPolymerTypes() == (PEPTIDE | CHEM));

    auto* clear_all =
        fixture.popup.findChild<QToolButton*>("filter_clear_all_link");
    BOOST_REQUIRE(clear_all != nullptr);
    clear_all->click();
    BOOST_REQUIRE(emitted.size() == 4);
    BOOST_CHECK(emitted.back() == MonomerDatabaseFilter());
    fixture.checkCombos("Any", false, "Any", false, false, "Any or None", false,
                        false);
    BOOST_TEST(fixture.filters_applied_lbl->text() == "No filters applied");

    // Clear All forgets remembered values, so restoring Peptide brings back
    // neither Terminal nor None
    fixture.setPolymerTypeChecked("Nucleic acid", false);
    fixture.setPolymerTypeChecked("Structure", false);
    fixture.checkCombos("Peptide", true, "Any", true, false, "Any or None",
                        true, false);
}

/**
 * The user can uncheck the last polymer type in the dropdown, which shows
 * "None" and disables the other combo boxes, so that switching from one
 * polymer type to another can uncheck the old type before checking the new one
 */
BOOST_AUTO_TEST_CASE(no_polymer_types)
{
    FilterPopupFixture fixture{MonomerDatabaseFilter(PEPTIDE)};
    auto* combo = fixture.polymer_type_combo;
    combo->showPopup();
    auto* view = combo->view();
    auto click_item = [&](int row) {
        QTest::mouseClick(
            view->viewport(), Qt::LeftButton, Qt::NoModifier,
            view->visualRect(combo->model()->index(row, 0)).center());
    };
    click_item(0);
    BOOST_CHECK(fixture.popup.getFilter().getPolymerTypes() ==
                MonomerDbPolymerType::NONE);
    BOOST_TEST(view->isVisible());
    fixture.checkCombos("None", true, "Any", false, false, "Any or None", false,
                        false);
    BOOST_TEST(fixture.filters_applied_lbl->text() == "1 filter applied");

    click_item(1);
    BOOST_CHECK(fixture.popup.getFilter().getPolymerTypes() == NUCLEIC_ACID);
    fixture.checkCombos("Nucleic acid", true, "Any", false, false,
                        "Any or None", true, false);
    combo->hidePopup();
}

/**
 * Values hidden by applicability are shown as the permissive values, and
 * restored when they apply again for the same polymer type, including after
 * passing through no polymer types. Switching to Nucleic acid forgets them,
 * whether Nucleic acid is checked before or after Peptide is unchecked.
 */
BOOST_AUTO_TEST_CASE(remembered_values)
{
    FilterPopupFixture fixture{MonomerDatabaseFilter(PEPTIDE)};
    fixture.choose(fixture.monomer_type_combo, "Backbone");
    fixture.choose(fixture.natural_analog_combo, "C");

    // checking and unchecking Nucleic acid in the open dropdown restores both
    auto* combo = fixture.polymer_type_combo;
    combo->showPopup();
    auto* view = combo->view();
    auto click_item = [&](int row) {
        QTest::mouseClick(
            view->viewport(), Qt::LeftButton, Qt::NoModifier,
            view->visualRect(combo->model()->index(row, 0)).center());
    };
    click_item(1);
    BOOST_TEST(view->isVisible());
    fixture.checkCombos("Peptide, Nucleic acid", true, "Any", false, false,
                        "Any or None", false, false);
    BOOST_TEST(fixture.filters_applied_lbl->text() == "1 filter applied");
    click_item(1);
    fixture.checkCombos("Peptide", true, "Backbone", true, true, "C", true,
                        true);
    BOOST_TEST(fixture.filters_applied_lbl->text() == "3 filters applied");
    combo->hidePopup();

    // unchecking and rechecking Peptide restores both
    fixture.setPolymerTypeChecked("Peptide", false);
    fixture.checkCombos("None", true, "Any", false, false, "Any or None", false,
                        false);
    fixture.setPolymerTypeChecked("Peptide", true);
    fixture.checkCombos("Peptide", true, "Backbone", true, true, "C", true,
                        true);

    // unchecking Peptide and then checking Nucleic acid forgets both
    fixture.setPolymerTypeChecked("Peptide", false);
    fixture.setPolymerTypeChecked("Nucleic acid", true);
    fixture.checkCombos("Nucleic acid", true, "Any", false, false,
                        "Any or None", true, false);
    fixture.setPolymerTypeChecked("Nucleic acid", false);
    fixture.setPolymerTypeChecked("Peptide", true);
    fixture.checkCombos("Peptide", true, "Any", true, false, "Any or None",
                        true, false);
    fixture.choose(fixture.monomer_type_combo, "Backbone");
    fixture.choose(fixture.natural_analog_combo, "C");

    // checking Nucleic acid and then unchecking Peptide forgets both
    fixture.setPolymerTypeChecked("Nucleic acid", true);
    fixture.setPolymerTypeChecked("Peptide", false);
    fixture.checkCombos("Nucleic acid", true, "Any", false, false,
                        "Any or None", true, false);
    fixture.setPolymerTypeChecked("Peptide", true);
    fixture.setPolymerTypeChecked("Nucleic acid", false);
    fixture.checkCombos("Peptide", true, "Any", true, false, "Any or None",
                        true, false);
}

/**
 * Hiding the popup closes any open dropdown
 */
BOOST_AUTO_TEST_CASE(hide_closes_dropdown)
{
    FilterPopupFixture fixture{MonomerDatabaseFilter(PEPTIDE)};
    fixture.polymer_type_combo->showPopup();
    BOOST_REQUIRE(fixture.polymer_type_combo->view()->isVisible());
    fixture.popup.hide();
    BOOST_TEST(!fixture.polymer_type_combo->view()->isVisible());
}

} // namespace sketcher
} // namespace schrodinger
