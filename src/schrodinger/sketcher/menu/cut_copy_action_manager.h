#pragma once

#include <QAction>
#include <QMenu>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/definitions.h"

enum class SceneSubset;

namespace schrodinger
{
namespace sketcher
{

class SketcherModel;

class SKETCHER_API CutCopyActionManager : public QWidget
{
    Q_OBJECT
  public:
    CutCopyActionManager(QWidget* parent);
    void setModel(SketcherModel* model);

    /**
     * @param whether copy actions should dynamically change based on selection
     */
    void setAlwaysCopyAll(bool always_copy_all);

    QAction* m_cut_action = nullptr;
    QAction* m_copy_action = nullptr;
    QMenu* m_copy_as_menu = nullptr;

  signals:
    void cutRequested(schrodinger::rdkit_extensions::Format format);
    void copyRequested(schrodinger::rdkit_extensions::Format format,
                       SceneSubset subset);

  private:
    /**
     * Updates action labels based on what subset is emitted on cut/copy
     */
    void updateActions();

    /**
     * @return subset to emit on cut/copy requests
     */
    SceneSubset getSubset();

    /**
     * Initializes the Copy As menu with both standard and reaction formats
     */
    void initCopyAsMenu();

    SketcherModel* m_sketcher_model = nullptr;
    bool m_always_copy_all = false;
};

} // namespace sketcher
} // namespace schrodinger
