#pragma once

#include <QToolButton>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A QToolButton subclass used for the AMINO/NUCLEIC toggle buttons in the
 * monomer tool widget. Applies the appropriate style sheet and a QProxyStyle
 * that prevents the label from shifting when the button is checked.
 */
class SKETCHER_API AminoOrNucleicToggleButton : public QToolButton
{
    Q_OBJECT
  public:
    AminoOrNucleicToggleButton(QWidget* parent = nullptr);
};

} // namespace sketcher
} // namespace schrodinger
