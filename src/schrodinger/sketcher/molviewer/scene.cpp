#include "schrodinger/sketcher/molviewer/scene.h"

#include <functional>
#include <unordered_set>

#include <GraphMol/ROMol.h>
#include <QApplication>
#include <QFont>
#include <QGraphicsSceneMouseEvent>
#include <QMimeData>
#include <QPainterPath>
#include <QString>
#include <QTransform>
#include <QUndoStack>
#include <QUrl>
#include <QtGlobal>

#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/qt_utils.h"

#define SETTER_AND_GETTER(settings_member, update_method, type, getter, \
                          setter, variable_name)                        \
    type Scene::getter() const                                          \
    {                                                                   \
        return settings_member.variable_name;                           \
    }                                                                   \
    void Scene::setter(type value)                                      \
    {                                                                   \
        settings_member.variable_name = value;                          \
        update_method();                                                \
    }
#define ATOM_SETTING(type, getter, setter, variable_name)                 \
    SETTER_AND_GETTER(m_atom_item_settings, updateAtomAndBondItems, type, \
                      getter, setter, variable_name)
#define BOND_SETTING(type, getter, setter, variable_name)                  \
    SETTER_AND_GETTER(m_bond_item_settings, updateBondItems, type, getter, \
                      setter, variable_name)

namespace schrodinger
{
namespace sketcher
{

NullSceneTool::NullSceneTool() : AbstractSceneTool(nullptr, nullptr)
{
}

Scene::Scene(MolModel* mol_model, SketcherModel* sketcher_model,
             QObject* parent) :
    QGraphicsScene(parent),
    m_mol_model(mol_model),
    m_sketcher_model(sketcher_model)
{
    m_selection_highlighting_item = new SelectionHighlightingItem();
    m_scene_tool = std::make_shared<NullSceneTool>();

    addItem(m_selection_highlighting_item);

    if (m_mol_model == nullptr || m_sketcher_model == nullptr) {
        throw std::runtime_error("Cannot construct without valid models");
    }
    connect(this, &Scene::changed, m_sketcher_model,
            &SketcherModel::interactiveItemsChanged);
    connect(this, &Scene::selectionChanged, m_sketcher_model,
            &SketcherModel::selectionChanged);

    connect(m_mol_model, &MolModel::moleculeChanged, this,
            &Scene::updateInteractiveItems);
    connect(m_mol_model, &MolModel::selectionChanged, this,
            &Scene::onMolModelSelectionChanged);

    connect(m_sketcher_model, &SketcherModel::valuesChanged, this,
            &Scene::onModelValuesChanged);
    connect(m_sketcher_model, &SketcherModel::interactiveItemsRequested, this,
            [this]() { return getInteractiveItems(); });
    connect(m_sketcher_model, &SketcherModel::selectionRequested, this,
            [this]() { return selectedItems(); });
    connect(m_sketcher_model, &SketcherModel::contextMenuObjectsRequested, this,
            [this]() { return m_context_menu_objects; });

    updateSceneTool();
}

Scene::~Scene()
{
    // load an empty scene tool to ensure that we unload the current scene tool.
    // Otherwise, both this class and the scene tool itself will try to destroy
    // any graphics items owned by the scene tool.
    std::shared_ptr<AbstractSceneTool> null_scene_tool =
        std::make_shared<NullSceneTool>();
    setSceneTool(null_scene_tool);

    // Avoid selection update calls when any selected item is destroyed
    disconnect(this, &Scene::selectionChanged, this,
               &Scene::updateSelectionHighlighting);
}

void Scene::updateInteractiveItems()
{
    clearAllInteractiveItems();
    const auto* mol = m_mol_model->getMol();

    m_atom_to_atom_item.clear();
    m_bond_to_bond_item.clear();

    // create atom items
    for (std::size_t i = 0; i < mol->getNumAtoms(); ++i) {
        const auto* atom = mol->getAtomWithIdx(i);
        const auto pos = mol->getConformer().getAtomPos(i);
        auto* atom_item = new AtomItem(atom, m_fonts, m_atom_item_settings);
        atom_item->setPos(to_scene_xy(pos));
        m_atom_to_atom_item[atom] = atom_item;
        addItem(atom_item);
    }

    // create bond items
    for (auto bond : mol->bonds()) {
        const auto* from_atom_item = m_atom_to_atom_item[bond->getBeginAtom()];
        const auto* to_atom_item = m_atom_to_atom_item[bond->getEndAtom()];
        auto* bond_item = new BondItem(bond, *from_atom_item, *to_atom_item,
                                       m_fonts, m_bond_item_settings);
        m_bond_to_bond_item[bond] = bond_item;
        addItem(bond_item);
    }
}

QList<QGraphicsItem*> Scene::getInteractiveItems() const
{
    // We expect all objects that inherit from AbstractGraphicsItem to
    // be interactive (atoms, bonds, etc.), and all objects that don't
    // to be purely graphical (selection highlighting paths, etc.)
    QList<QGraphicsItem*> interactive_items;
    for (auto item : items()) {
        if (dynamic_cast<AbstractGraphicsItem*>(item) != nullptr) {
            interactive_items.append(item);
        }
    }
    return interactive_items;
}

void Scene::clearAllInteractiveItems()
{
    // remove all interactive items and reset the rdkit molecule; this will
    // preserve items include selection paths, highlighting items, and
    // potentially the watermark managed by the SketcherWidget
    for (auto item : getInteractiveItems()) {
        removeItem(item);
        delete item;
    }
}

qreal Scene::fontSize() const
{
    return m_fonts.size();
}

void Scene::setFontSize(qreal size)
{
    m_fonts.setSize(size);
    updateAtomAndBondItems();
}

ATOM_SETTING(CarbonLabels, carbonsLabeled, setCarbonsLabeled, m_carbon_labels)
ATOM_SETTING(bool, valenceErrorsShown, setValenceErrorsShown,
             m_valence_errors_shown)
BOND_SETTING(qreal, bondWidth, setBondWidth, m_bond_width)
BOND_SETTING(qreal, doubleBondSpacing, setDoubleBondSpacing,
             m_double_bond_spacing)

void Scene::updateAtomAndBondItems()
{
    // we need to update all of the atom items before we update any bond items,
    // since bond items pull information from their associated atom items
    for (auto item : items()) {
        if (auto atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            atom_item->updateCachedData();
        }
    }
    updateBondItems();
}

void Scene::updateBondItems()
{
    for (auto item : items()) {
        if (auto bond_item = qgraphicsitem_cast<BondItem*>(item)) {
            bond_item->updateCachedData();
        }
    }
}

void Scene::updateSelectionHighlighting()
{
    m_selection_highlighting_item->highlightItems(selectedItems());
}

void Scene::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    m_mouse_down_screen_pos = event->screenPos();
    if (event->button() == Qt::LeftButton) {
        m_scene_tool->onMousePress(event);
    }
}

void Scene::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // TODO: don't pass events to the scene tool if a context menu is visible
    // TODO: SKETCH-1890: translate on right click and rotate on middle click
    m_scene_tool->onMouseMove(event);
    if (event->buttons() & Qt::LeftButton) {
        // check to see whether the mouse has moved far enough to start a
        // drag
        if (!m_drag_started) {
            int drag_dist = (m_mouse_down_screen_pos - event->screenPos())
                                .manhattanLength();
            if (drag_dist >= QApplication::startDragDistance()) {
                m_drag_started = true;
                m_scene_tool->onDragStart(event);
            }
        }
        if (m_drag_started) {
            m_scene_tool->onDragMove(event);
        }
    }
}

void Scene::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    if (event->button() == Qt::LeftButton) {
        if (m_drag_started) {
            m_scene_tool->onDragRelease(event);
        } else {
            m_scene_tool->onMouseClick(event);
        }
    }

    m_mouse_down_screen_pos = QPointF();
    m_drag_started = false;
}

void Scene::selectGraphicsItems(const QList<QGraphicsItem*>& items,
                                const SelectMode select_mode)
{
    std::unordered_set<const RDKit::Atom*> atoms;
    std::unordered_set<const RDKit::Bond*> bonds;
    for (auto cur_item : items) {
        if (auto atom_item = qgraphicsitem_cast<AtomItem*>(cur_item)) {
            atoms.insert(atom_item->getAtom());
        } else if (auto bond_item = qgraphicsitem_cast<BondItem*>(cur_item)) {
            bonds.insert(bond_item->getBond());
        }
    }
    m_mol_model->select(atoms, bonds, select_mode);
}

void Scene::onMolModelSelectionChanged()
{
    // block the per-item signals for performance reasons
    Scene::SelectionChangeSignalBlocker signal_blocker(this);
    clearSelection();
    for (auto* atom : m_mol_model->getSelectedAtoms()) {
        m_atom_to_atom_item[atom]->setSelected(true);
    }
    for (auto* bond : m_mol_model->getSelectedBonds()) {
        m_bond_to_bond_item[bond]->setSelected(true);
    }
    updateSelectionHighlighting();
    // signal_blocker emits selectionChanged on destruction
}

void Scene::onModelValuesChanged(const std::unordered_set<ModelKey>& keys)
{
    for (auto key : keys) {
        switch (key) {
            case ModelKey::SELECTION_TOOL:
            case ModelKey::DRAW_TOOL:
            case ModelKey::ATOM_TOOL:
            case ModelKey::BOND_TOOL:
            case ModelKey::CHARGE_TOOL:
            case ModelKey::RING_TOOL:
            case ModelKey::ENUMERATION_TOOL:
            case ModelKey::ELEMENT:
            case ModelKey::ATOM_QUERY:
            case ModelKey::RGROUP_NUMBER:
                updateSceneTool();
                break;
            default:
                break;
        }
    }
}

void Scene::updateSceneTool()
{
    std::shared_ptr<AbstractSceneTool> new_scene_tool;
    auto draw_tool = m_sketcher_model->getDrawTool();
    if (draw_tool == DrawTool::SELECT) {
        auto select_tool = m_sketcher_model->getSelectionTool();
        new_scene_tool = get_select_scene_tool(select_tool, this, m_mol_model);
    } else {
        // tool not yet implemented
        new_scene_tool = std::make_shared<NullSceneTool>();
    }
    setSceneTool(new_scene_tool);
}

void Scene::setSceneTool(std::shared_ptr<AbstractSceneTool> new_scene_tool)
{
    // remove any graphics items associated with the old scene tool
    for (auto* item : m_scene_tool->getGraphicsItems()) {
        removeItem(item);
    }
    m_scene_tool = new_scene_tool;
    // add graphics items from the new scene tool
    for (auto* item : m_scene_tool->getGraphicsItems()) {
        addItem(item);
    }
}

Scene::SelectionChangeSignalBlocker::SelectionChangeSignalBlocker(
    Scene* scene) :
    m_scene(scene),
    m_original_value(m_scene->signalsBlocked())
{
    m_scene->blockSignals(true);
}
Scene::SelectionChangeSignalBlocker::~SelectionChangeSignalBlocker()
{
    m_scene->blockSignals(m_original_value);
    m_scene->selectionChanged();
}

QWidget* Scene::window() const
{
    auto parent_widget = dynamic_cast<QWidget*>(parent());
    return parent_widget ? parent_widget->window() : nullptr;
}

void Scene::dropEvent(QGraphicsSceneDragDropEvent* event)
{
    auto mime_data = event->mimeData();
    for (auto url : mime_data->urls()) {
        if (url.isLocalFile()) { // File drop event
            auto file_path = url.toLocalFile();
            auto contents = get_file_text(file_path.toStdString());
            auto format = get_file_format(file_path);
            emit importTextRequested(contents, format);
        }
    }
    event->acceptProposedAction();
}

void Scene::dragEnterEvent(QGraphicsSceneDragDropEvent* event)
{
    event->acceptProposedAction();
}

void Scene::dragMoveEvent(QGraphicsSceneDragDropEvent* event)
{
    event->acceptProposedAction();
}

void Scene::dragLeaveEvent(QGraphicsSceneDragDropEvent* event)
{
    event->acceptProposedAction();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/scene.moc"
