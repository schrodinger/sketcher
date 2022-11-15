/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <QGraphicsView>

#include "schrodinger/sketcher/definitions.h"

class QGraphicsScene;
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
};

} // namespace sketcher
} // namespace schrodinger
