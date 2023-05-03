#pragma once
#include <QMenu>

#include "schrodinger/sketcher/Bond.h"
#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API ModifyBondsMenu : public QMenu
{
    Q_OBJECT
  public:
    ModifyBondsMenu(QWidget* parent = nullptr);
    void setFlipEnabled(bool b);

  signals:
    void changeTypeRequested(sketcherBond::BondType type);
    void flipRequested();

  private:
    QMenu* createOtherTypeMenu();
    QMenu* createQueryMenu();
    QAction* m_flip_action;
};

class BondContextMenu : public ModifyBondsMenu
{
    Q_OBJECT
  public:
    BondContextMenu(QWidget* parent = nullptr);

  signals:
    void deleteRequested();
};

} // namespace sketcher
} // namespace schrodinger
