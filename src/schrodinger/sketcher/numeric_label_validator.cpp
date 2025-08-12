#include "schrodinger/sketcher/numeric_label_validator.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * Validator for the "numeric label" line edit of the bracket subgroup dialog
 */

NumericLabelValidator::NumericLabelValidator(QObject* parent) :
    QIntValidator(parent)
{
    setRange(0, 99999);
}

QValidator::State NumericLabelValidator::validate(QString& input,
                                                  int& pos) const
{
    if (QIntValidator::validate(input, pos) == QValidator::Acceptable) {
        return QValidator::Acceptable;
    }

    input = input.trimmed();

    if (input == "n") {
        return QValidator::Acceptable;
    }

    QString range_str{"-"};
    auto substrings = input.split(range_str);
    auto validity = substrings.size() <= 2 ? QValidator::Acceptable
                                           : QValidator::Intermediate;
    int subpos = 0;
    for (auto substring : substrings) {
        auto sub_validity = QIntValidator::validate(substring, subpos);
        if (substring.isEmpty()) {
            validity = QValidator::Intermediate;
        }
        if (sub_validity == QValidator::Invalid) {
            return sub_validity;
        } else if (sub_validity == QValidator::Intermediate) {
            validity = sub_validity;
        }
    }
    return validity;
}

void NumericLabelValidator::fixup(QString& input) const
{
}

} // namespace sketcher
} // namespace schrodinger

#include "numeric_label_validator.moc"
