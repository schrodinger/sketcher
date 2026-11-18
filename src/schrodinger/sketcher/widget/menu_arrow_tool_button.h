#pragma once

#include <QToolButton>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

/**
 * A tool button that shows the platform style's standard menu arrow while it's
 * hovered or checked, to show that clicking it opens a menu or popup. Its icon
 * and text are centered together, ignoring the arrow.
 */
class SKETCHER_API MenuArrowToolButton : public QToolButton
{
    Q_OBJECT

  public:
    MenuArrowToolButton(QWidget* parent = nullptr);

  protected:
    /**
     * Overridden QToolButton method: Paint centered contents with a standard
     * menu arrow on hover or check.
     */
    void paintEvent(QPaintEvent* event) override;
};

} // namespace sketcher
} // namespace schrodinger
