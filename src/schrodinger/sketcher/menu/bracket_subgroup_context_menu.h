#pragma once

#include <QMenu>

#include "schrodinger/sketcher/definitions.h"

namespace schrodinger
{
namespace sketcher
{

class SKETCHER_API BracketSubgroupContextMenu : public QMenu
{
    Q_OBJECT

  public:
    BracketSubgroupContextMenu(QWidget* parent = nullptr);

  signals:
    void bracketSubgroupDialogRequested() const;
    void deleteRequested() const;
};

} // namespace sketcher
} // namespace schrodinger
