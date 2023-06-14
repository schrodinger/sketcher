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
#include <QWidget>
#include <QtGlobal>
#include <QScreen>

#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/tool/select_erase_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_atom_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_bond_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_r_group_scene_tool.h"
#include "schrodinger/sketcher/tool/rotate_scene_tool.h"
#include "schrodinger/sketcher/tool/translate_scene_tool.h"

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
             QWidget* parent) :
    QGraphicsScene(parent),
    m_mol_model(mol_model),
    m_sketcher_model(sketcher_model)
{
    m_selection_highlighting_item = new SelectionHighlightingItem();
    m_left_button_scene_tool = std::make_shared<NullSceneTool>();
    m_middle_button_scene_tool =
        std::make_shared<RotateSceneTool>(this, m_mol_model);
    m_right_button_scene_tool =
        std::make_shared<TranslateSceneTool>(this, m_mol_model);

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

    m_background_context_menu =
        new BackgroundContextMenu(m_sketcher_model, parent);
    connectContextMenu(*m_background_context_menu);

    updateSceneTool();
}

void Scene::connectContextMenu(const BackgroundContextMenu& menu)
{
    /*
    connect(&menu, &BackgroundContextMenu::saveImageRequested, this,
            &sketcherScene::onSaveImageRequested);
    connect(&menu, &BackgroundContextMenu::exportToFileRequested, this,
            &sketcherScene::exportToFileRequested);
    connect(&menu, &BackgroundContextMenu::undoRequested, this,
            &sketcherScene::undo);
    connect(&menu, &BackgroundContextMenu::redoRequested, this,
            &sketcherScene::redo);
            */
    connect(&menu, &BackgroundContextMenu::flipHorizontalRequested, m_mol_model,
            &MolModel::flipAllHorizontal);
    connect(&menu, &BackgroundContextMenu::flipVerticalRequested, m_mol_model,
            &MolModel::flipAllVertical);
    /*
connect(&menu, &BackgroundContextMenu::selectAllRequested, this,
    &sketcherScene::selectAll);
connect(&menu, &BackgroundContextMenu::copyRequested, this,
    &sketcherScene::copyRequested);
connect(&menu, &BackgroundContextMenu::pasteRequested, this,
    [this]() { emit pasteRequested(true); });
connect(&menu, &BackgroundContextMenu::clearRequested, this,
    &sketcherScene::clearStructure);
connect(&menu, &BackgroundContextMenu::aboutToHide, this,
    &sketcherScene::onContextMenuHidden);
    */
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
    auto scene_tool = getSceneTool(event);
    scene_tool->onMousePress(event);
    event->accept();
}

void Scene::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // TODO: don't pass events to the scene tool if a context menu is
    // visible

    auto scene_tool = getSceneTool(event);
    if (event->buttons() != Qt::NoButton && !m_drag_started) {
        // check to see whether the mouse has moved far enough to start a drag
        int drag_dist =
            (m_mouse_down_screen_pos - event->screenPos()).manhattanLength();
        if (drag_dist >= QApplication::startDragDistance()) {
            m_drag_started = true;
            scene_tool->onDragStart(event);
        }
    }
    scene_tool->onMouseMove(event);
    /*make sure  m_left_button_scene_tool->onMouseMove gets always called,
     * regardless of the button which is actually being used, to update
     * predictive highlighting */
    if (scene_tool != m_left_button_scene_tool) {
        m_left_button_scene_tool->onMouseMove(event);
    }
    if (m_drag_started) {
        scene_tool->onDragMove(event);
    }

    event->accept();
}

void Scene::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    auto scene_tool = getSceneTool(event);
    if (m_drag_started) {
        scene_tool->onDragRelease(event);
    } else {
        scene_tool->onMouseClick(event);
    }
    scene_tool->onMouseRelease(event);
    m_mouse_down_screen_pos = QPointF();
    m_drag_started = false;
    event->accept();
}

void Scene::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event)
{
    // TODO: select clicked molecule
}

void Scene::showContextMenu(QGraphicsSceneMouseEvent* event)
{
    // Collect the set of items that the context menu should interact with
    // based on the position of the cursor
    std::unordered_set<QGraphicsItem*> context_menu_objects;
    auto top_item = getTopInteractiveItemAt(event->pos());
    if (top_item != nullptr) {
        if (top_item->isSelected()) {
            // If the cursor is over a selected object, the context menu
            // should be concerned with all selected objects
            for (auto item : selectedItems()) {
                context_menu_objects.insert(item);
            }
        } else {
            context_menu_objects.insert(top_item);
        }
    }
    // TODO: call different menus based on context menu objects. For now,
    // just call m_background_context_menu
    QMenu* menu = nullptr;
    menu = m_background_context_menu;

    menu->move(event->screenPos());
    auto screen_rect = QApplication::screenAt(QCursor::pos())->geometry();
    auto menu_rectangle = menu->geometry();

    // Make sure the menu is not off the screen (or on a different screen)
    if (menu_rectangle.left() < screen_rect.left()) {
        menu_rectangle.moveLeft(screen_rect.left());
    }
    if (menu_rectangle.top() < screen_rect.top()) {
        menu_rectangle.moveTop(screen_rect.top());
    }
    if (menu_rectangle.right() > screen_rect.right()) {
        menu_rectangle.moveRight(screen_rect.right());
    }
    if (menu_rectangle.bottom() > screen_rect.bottom()) {
        menu_rectangle.moveBottom(screen_rect.bottom());
    }
    menu->move(menu_rectangle.topLeft());
    menu->show();
}

AbstractGraphicsItem* Scene::getTopInteractiveItemAt(const QPointF& pos) const
{
    for (auto item : items(pos)) {
        if (auto interactive_item = dynamic_cast<AbstractGraphicsItem*>(item)) {
            return interactive_item;
        }
    }
    return nullptr;
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

std::shared_ptr<AbstractSceneTool>
Scene::getSceneTool(QGraphicsSceneMouseEvent* const event)
{
    // the button for a release event is stored in event->button(), while the
    // one for a drag event is stored in event->buttons()
    if ((event->buttons() & Qt::MiddleButton) ||
        (event->button() == Qt::MiddleButton)) {
        return m_middle_button_scene_tool;
    } else if ((event->buttons() & Qt::RightButton) ||
               (event->button() == Qt::RightButton)) {
        return m_right_button_scene_tool;
    } else {
        return m_left_button_scene_tool;
    }
}

void Scene::updateSceneTool()
{
    std::shared_ptr<AbstractSceneTool> new_scene_tool = getNewSceneTool();
    setSceneTool(new_scene_tool);
}

std::shared_ptr<AbstractSceneTool> Scene::getNewSceneTool()
{
    auto draw_tool = m_sketcher_model->getDrawTool();
    if (draw_tool == DrawTool::SELECT) {
        auto select_tool = m_sketcher_model->getSelectionTool();
        return get_select_scene_tool(select_tool, this, m_mol_model);
    } else if (draw_tool == DrawTool::ERASE) {
        return std::make_shared<EraseSceneTool>(this, m_mol_model);
    } else if (draw_tool == DrawTool::ATOM) {
        auto atom_tool = m_sketcher_model->getAtomTool();
        if (atom_tool == AtomTool::ELEMENT) {
            auto element = m_sketcher_model->getElement();
            return std::make_shared<DrawElementSceneTool>(element, this,
                                                          m_mol_model);
        } else { // atom_tool == AtomTool::AtomQuery
            auto atom_query = m_sketcher_model->getAtomQuery();
            return std::make_shared<DrawAtomQuerySceneTool>(atom_query, this,
                                                            m_mol_model);
        }
    } else if (draw_tool == DrawTool::BOND) {
        auto bond_tool = m_sketcher_model->getBondTool();
        switch (bond_tool) {
            case BondTool::SINGLE:
            case BondTool::DOUBLE:
            case BondTool::TRIPLE:
            case BondTool::COORDINATE:
            case BondTool::ZERO:
            case BondTool::SINGLE_UP:
            case BondTool::SINGLE_DOWN:
            case BondTool::AROMATIC:
            case BondTool::SINGLE_EITHER:
            case BondTool::DOUBLE_EITHER:
                return std::make_shared<DrawBondSceneTool>(bond_tool, this,
                                                           m_mol_model);
                break;
            case BondTool::SINGLE_OR_DOUBLE:
            case BondTool::SINGLE_OR_AROMATIC:
            case BondTool::DOUBLE_OR_AROMATIC:
            case BondTool::ANY:
                return std::make_shared<DrawBondQuerySceneTool>(bond_tool, this,
                                                                m_mol_model);
                break;
            case BondTool::ATOM_CHAIN:
                // TODO
                // return std::make_shared<DrawChainSceneTool>(this,
                // m_mol_model);
                break;
        }
    } else if (draw_tool == DrawTool::ENUMERATION) {
        auto enumeration_tool = m_sketcher_model->getEnumerationTool();
        if (enumeration_tool == EnumerationTool::NEW_RGROUP) {
            return std::make_shared<DrawIncrementingRGroupSceneTool>(
                this, m_mol_model);
        } else if (enumeration_tool == EnumerationTool::EXISTING_RGROUP) {
            auto r_group_num =
                m_sketcher_model->getValueInt(ModelKey::RGROUP_NUMBER);
            return std::make_shared<DrawRGroupSceneTool>(r_group_num, this,
                                                         m_mol_model);
        } else if (enumeration_tool == EnumerationTool::ATTACHMENT_POINT) {
            // TODO
            // return std::make_shared<DrawAttachmentPointSceneTool>(this,
            // m_mol_model);
        }
    }
    // tool not yet implemented
    return std::make_shared<NullSceneTool>();
}

void Scene::setSceneTool(std::shared_ptr<AbstractSceneTool> new_scene_tool)
{
    // remove any graphics items associated with the old scene tool
    for (auto* item : m_left_button_scene_tool->getGraphicsItems()) {
        removeItem(item);
    }
    m_left_button_scene_tool = new_scene_tool;
    // add graphics items from the new scene tool
    for (auto* item : m_left_button_scene_tool->getGraphicsItems()) {
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
