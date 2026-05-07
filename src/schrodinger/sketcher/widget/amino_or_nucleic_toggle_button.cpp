#include "schrodinger/sketcher/widget/amino_or_nucleic_toggle_button.h"

#include <QProxyStyle>

#include "schrodinger/sketcher/sketcher_css_style.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A QProxyStyle that prevents the button's label from shifting when the button
 * is checked. Normally the label is lowered so that the button looks like it's
 * been pressed in.
 */
class NoShiftProxyStyle : public QProxyStyle
{
  public:
    using QProxyStyle::QProxyStyle;

    int pixelMetric(PixelMetric metric, const QStyleOption* option = nullptr,
                    const QWidget* widget = nullptr) const override
    {
        if (metric == PM_ButtonShiftHorizontal ||
            metric == PM_ButtonShiftVertical) {
            return 0;
        }
        return QProxyStyle::pixelMetric(metric, option, widget);
    }
};

AminoOrNucleicToggleButton::AminoOrNucleicToggleButton(QWidget* parent) :
    QToolButton(parent)
{
    setStyleSheet(AMINO_OR_NUCLEIC_TOGGLE_STYLE);
    auto* proxy_style = new NoShiftProxyStyle();
    proxy_style->setParent(this);
    setStyle(proxy_style);
}

} // namespace sketcher
} // namespace schrodinger
