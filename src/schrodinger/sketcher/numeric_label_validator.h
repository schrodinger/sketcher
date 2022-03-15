#pragma once

#include "schrodinger/sketcher/definitions.h"
#include <QIntValidator>

namespace schrodinger
{
namespace sketcher
{

/**
 * Validator for the "numeric label" line edit of the bracket subgroup dialog
 */
class SKETCHER_API NumericLabelValidator : public QIntValidator
{
    Q_OBJECT

  public:
    NumericLabelValidator(QObject* parent = nullptr);
    QValidator::State validate(QString& input, int& pos) const override;
    void fixup(QString& input) const override;
};

} // namespace sketcher
} // namespace schrodinger
