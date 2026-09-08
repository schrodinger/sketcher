#define BOOST_TEST_MODULE test_monomer_tool_widget
#include <boost/test/unit_test.hpp>

#include <QAbstractButton>

#include "../test_common.h"
#include "schrodinger/sketcher/widget/monomer_tool_widget.h"

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

namespace schrodinger
{
namespace sketcher
{

/**
 * Verify that unknown monomer buttons use the unknown monomer styling. If these
 * buttons are converted to ModularToolButtons, then the ModularToolButton
 * styling will inadvertently overwrite the UNKNOWN_MONOMER_STYLE.
 */
BOOST_AUTO_TEST_CASE(unknown_monomer_button_styles)
{
    MonomerToolWidget widget;

    for (const auto* button_name : {"unk_btn", "na_n_btn"}) {
        const auto* button = widget.findChild<QAbstractButton*>(button_name);
        BOOST_REQUIRE(button != nullptr);
        const auto style_sheet = button->styleSheet();
        BOOST_TEST_CONTEXT(button->objectName().toStdString())
        {
            BOOST_TEST(style_sheet.contains(QStringLiteral("italic")));
            BOOST_TEST(style_sheet.contains(QStringLiteral("#606060")));
        }
    }
}

} // namespace sketcher
} // namespace schrodinger
