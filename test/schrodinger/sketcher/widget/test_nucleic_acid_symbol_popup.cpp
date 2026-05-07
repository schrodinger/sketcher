#define BOOST_TEST_MODULE Test_Sketcher
#include <algorithm>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/widget/monomer_tool_widget.h"
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

namespace
{
// Subclass exposes the protected popup map so we can inspect what the
// MonomerToolWidget constructor wired up for each NA button.
class MonomerToolWidgetForTest : public MonomerToolWidget
{
  public:
    using MonomerToolWidget::m_nucleic_acid_symbol_popups;
};

static std::vector<std::string>
get_analog_symbols(const NucleicAcidSymbolPopup& popup)
{
    std::vector<std::string> symbols;
    for (const auto& packet : popup.getButtonPackets()) {
        if (packet.enum_int == 0) {
            continue; // skip the standard at id 0
        }
        symbols.push_back(popup.getSymbolForId(packet.enum_int).toStdString());
    }
    return symbols;
}

static NucleicAcidSymbolPopup*
find_popup_for_standard(const MonomerToolWidgetForTest& widget,
                        const std::string& standard_symbol)
{
    for (const auto& [btn, popup] : widget.m_nucleic_acid_symbol_popups) {
        if (popup->getSymbolForId(0).toStdString() == standard_symbol) {
            return popup;
        }
    }
    return nullptr;
}
} // namespace

/**
 * dR has NATURAL_ANALOG=R, so the analog map has no entry under "dR". The
 * widget should pull R's analog list (and R itself) into dR's popup so it
 * mirrors R's popup, which has long shown dR.
 */
BOOST_AUTO_TEST_CASE(test_dr_popup_inherits_r_analogs)
{
    MonomerToolWidgetForTest widget;
    auto* dr_popup = find_popup_for_standard(widget, "dR");
    auto* r_popup = find_popup_for_standard(widget, "R");
    BOOST_REQUIRE(dr_popup != nullptr);
    BOOST_REQUIRE(r_popup != nullptr);

    auto dr_analogs = get_analog_symbols(*dr_popup);
    auto r_analogs = get_analog_symbols(*r_popup);

    auto contains = [](const std::vector<std::string>& v,
                       const std::string& s) {
        return std::find(v.begin(), v.end(), s) != v.end();
    };
    // R appears as a clickable analog in dR's popup (the parent of dR).
    BOOST_TEST(contains(dr_analogs, "R"));
    // dR is the standard at id 0, not duplicated as an analog of itself.
    BOOST_TEST(!contains(dr_analogs, "dR"));
    // R's popup still shows dR (existing symmetric behavior).
    BOOST_TEST(contains(r_analogs, "dR"));
}

} // namespace sketcher
} // namespace schrodinger
