
#define BOOST_TEST_MODULE Test_Sketcher
#include <boost/test/unit_test.hpp>

#include "schrodinger/sketcher/numeric_label_validator.h"
#include "test_common.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Verify that the `validate()` method works as expected.
 */
BOOST_AUTO_TEST_CASE(test_validate)
{
    NumericLabelValidator validator;
    int pos = 0;

    for (QString input : {"n", "100", "0", "99-4", "01234"}) {
        BOOST_TEST(validator.validate(input, pos) == QValidator::Acceptable,
                   input.toStdString());
    }

    for (QString input : {"1-2-3", "23--", "-1", "-"}) {
        BOOST_TEST(validator.validate(input, pos) == QValidator::Intermediate,
                   input.toStdString());
    }

    for (QString input : {"34abc", "nn", "100-200-x-400"}) {
        BOOST_TEST(validator.validate(input, pos) == QValidator::Invalid,
                   input.toStdString());
    }
}

} // namespace sketcher
} // namespace schrodinger
