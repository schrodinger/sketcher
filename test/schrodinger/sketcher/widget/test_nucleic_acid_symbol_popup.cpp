#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/widget/nucleic_acid_symbol_popup.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

static rdkit_extensions::MonomerInfo make_analog(const std::string& symbol,
                                                 const std::string& name)
{
    rdkit_extensions::MonomerInfo info;
    info.symbol = symbol;
    info.name = name;
    return info;
}

/**
 * Construct with no analogs: only the standard base button is present at ID 0.
 */
BOOST_AUTO_TEST_CASE(test_standard_only)
{
    NucleicAcidSymbolPopup popup("C", "Cytosine", {});
    auto packets = popup.getButtonPackets();
    BOOST_TEST(packets.size() == 1);
    BOOST_TEST(popup.getSymbolForId(0).toStdString() == "C");
    BOOST_TEST(popup.getSymbolForId(1).isEmpty());
}

/**
 * Construct with two analogs: standard at ID 0, analogs at IDs 1 and 2.
 */
BOOST_AUTO_TEST_CASE(test_with_analogs)
{
    std::vector<rdkit_extensions::MonomerInfo> analogs = {
        make_analog("5mC", "5-methylcytosine"),
        make_analog("hC", "hydroxymethylcytosine"),
    };
    NucleicAcidSymbolPopup popup("C", "Cytosine", analogs);
    auto packets = popup.getButtonPackets();
    BOOST_TEST(packets.size() == 3);
    BOOST_TEST(popup.getSymbolForId(0).toStdString() == "C");
    BOOST_TEST(popup.getSymbolForId(1).toStdString() == "5mC");
    BOOST_TEST(popup.getSymbolForId(2).toStdString() == "hC");
    BOOST_TEST(popup.getSymbolForId(99).isEmpty());
}

/**
 * Each popup button carries a tooltip set to the human-readable monomer name,
 * which is what feeds the toolbar tooltip via ModularPopup::getToolTip().
 */
BOOST_AUTO_TEST_CASE(test_button_tooltips)
{
    std::vector<rdkit_extensions::MonomerInfo> analogs = {
        make_analog("5mC", "5-methylcytosine"),
    };
    NucleicAcidSymbolPopup popup("C", "Cytosine", analogs);
    auto packets = popup.getButtonPackets();
    BOOST_REQUIRE(packets.size() == 2);
    BOOST_TEST(packets[0].button->toolTip().toStdString() == "Cytosine");
    BOOST_TEST(packets[1].button->toolTip().toStdString() ==
               "5-methylcytosine");
}

} // namespace sketcher
} // namespace schrodinger
