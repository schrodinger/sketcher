#pragma once

#include <QGraphicsView>

#include "schrodinger/sketcher/definitions.h"

class QGraphicsScene;
class QResizeEvent;
class QWidget;

namespace schrodinger
{
namespace sketcher
{

/**
 * A Qt graphics view for displaying molecules in a molviewer Scene.
 */
class SKETCHER_API View : public QGraphicsView
{
    Q_OBJECT
  public:
    View(QGraphicsScene* scene, QWidget* parent = nullptr);
    View(QWidget* parent = nullptr);

  public slots:
    /**
     * fit the scene to the view so that the entire scene is visible, but avoid
     * zooming in too much for small molecules
     */
    void onFitToScreenRequested();

  signals:
    void resized();

  protected:
    // Override the QGraphicsView method so we can call enlargeSceneIfNeeded
    void resizeEvent(QResizeEvent* event) override;

    /**
     * Make sure that the scene is at least as large as the view.  If the scene
     * is smaller than the view, then the view will center the scene.  Then, if
     * a new item gets added that causes the scene to grow, the view will
     * re-center, which causes the molecule to jump around.
     */
    void enlargeSceneIfNeeded();

    /**
     * Make sure that the scene has enough space around the items so that it can
     * be centered
     */
    void adjustSceneAroundItems();
};

} // namespace sketcher
} // namespace schrodinger
