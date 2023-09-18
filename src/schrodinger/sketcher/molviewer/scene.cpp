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
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/non_molecular_item.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/tool/arrow_plus_scene_tool.h"
#include "schrodinger/sketcher/tool/atom_mapping_scene_tool.h"
#include "schrodinger/sketcher/tool/attachment_point_scene_tool.h"
#include "schrodinger/sketcher/tool/select_erase_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_atom_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_bond_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_chain_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_fragment_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_r_group_scene_tool.h"
#include "schrodinger/sketcher/tool/edit_charge_scene_tool.h"

#include "schrodinger/sketcher/molviewer/rotation_item.h"
#include "schrodinger/sketcher/tool/rotate_scene_tool.h"
#include "schrodinger/sketcher/tool/translate_scene_tool.h"
#include "schrodinger/sketcher/tool/move_rotate_scene_tool.h"

namespace schrodinger
{
namespace sketcher
{

namespace
{

/**
 * @return whether a graphics item is one of the types specified by the type
 * flag.
 */
bool item_matches_type_flag(QGraphicsItem* item,
                            InteractiveItemFlagType type_flag)
{
    auto type = item->type();
    if (type == AtomItem::Type) {
        if ((type_flag & InteractiveItemFlag::ATOM) ==
            InteractiveItemFlag::ATOM) {
            // we want all atoms, regardless of whether or not they're
            // attachment points
            return true;
        } else if (type_flag & InteractiveItemFlag::ATOM) {
            auto* atom = static_cast<AtomItem*>(item)->getAtom();
            if (is_r_group(atom)) {
                return type_flag & InteractiveItemFlag::R_GROUP;
            } else if (is_attachment_point(atom)) {
                return type_flag & InteractiveItemFlag::ATTACHMENT_POINT;
            } else {
                return type_flag & InteractiveItemFlag::ATOM_NOT_R_NOT_AP;
            }
        }
    } else if (type == BondItem::Type) {
        if ((type_flag & InteractiveItemFlag::BOND) ==
            InteractiveItemFlag::BOND) {
            // we want all bonds, regardless of whether or not they're
            // attachment point bonds
            return true;
        } else if (type_flag & InteractiveItemFlag::BOND) {
            auto* bond = static_cast<BondItem*>(item)->getBond();
            return is_attachment_point_bond(bond) ==
                   static_cast<bool>(
                       type_flag & InteractiveItemFlag::ATTACHMENT_POINT_BOND);
        }
    } else if (type == NonMolecularItem::Type) {
        return type_flag & InteractiveItemFlag::NON_MOLECULAR;
    }
    return false;
}

} // namespace

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
            &Scene::updateMolecularItems);
    connect(m_mol_model, &MolModel::nonMolecularObjectsChanged, this,
            &Scene::updateNonMolecularItems);

    connect(m_mol_model, &MolModel::coordinatesChanged, this,
            &Scene::moveInteractiveItems);

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

void Scene::moveInteractiveItems()
{

    update_conf_for_mol_graphics_items(
        getInteractiveItems(InteractiveItemFlag::ATOM),
        getInteractiveItems(InteractiveItemFlag::BOND), *m_mol_model->getMol());
    for (auto* non_mol_obj : m_mol_model->getNonMolecularObjects()) {
        auto* non_mol_item =
            m_non_molecular_to_non_molecular_item.at(non_mol_obj);
        non_mol_item->setPos(to_scene_xy(non_mol_obj->getCoords()));
        non_mol_item->updateCachedData();
    }
    updateSelectionHighlighting();
    m_left_button_scene_tool->onStructureUpdated();
}

void Scene::updateMolecularItems()
{
    clearInteractiveItems(InteractiveItemFlag::MOLECULAR);
    const auto* mol = m_mol_model->getMol();
    std::vector<QGraphicsItem*> all_items;
    std::tie(all_items, m_atom_to_atom_item, m_bond_to_bond_item) =
        create_graphics_items_for_mol(mol, m_fonts, m_atom_item_settings,
                                      m_bond_item_settings);
    for (auto* item : all_items) {
        addItem(item);
        m_interactive_items.insert(item);
    }
    updateItemSelection();
    updateSelectionHighlighting();
    m_left_button_scene_tool->onStructureUpdated();
}

void Scene::updateNonMolecularItems()
{
    clearInteractiveItems(InteractiveItemFlag::NON_MOLECULAR);
    QColor color = m_atom_item_settings.getAtomColor(-1);
    for (const auto* model_obj : m_mol_model->getNonMolecularObjects()) {
        auto* item = new NonMolecularItem(model_obj, color);
        auto rd_coords = model_obj->getCoords();
        auto qpos = to_scene_xy(rd_coords);
        item->setPos(qpos);
        m_non_molecular_to_non_molecular_item[model_obj] = item;
        addItem(item);
        m_interactive_items.insert(item);
    }
    updateItemSelection();
    updateSelectionHighlighting();
    m_left_button_scene_tool->onStructureUpdated();
}

void Scene::updateItemSelection()
{
    clearSelection();
    for (auto* atom : m_mol_model->getSelectedAtoms()) {
        m_atom_to_atom_item[atom]->setSelected(true);
    }
    for (auto* bond : m_mol_model->getSelectedBonds()) {
        m_bond_to_bond_item[bond]->setSelected(true);
    }
    for (auto* non_molecular_object :
         m_mol_model->getSelectedNonMolecularObjects()) {
        m_non_molecular_to_non_molecular_item[non_molecular_object]
            ->setSelected(true);
    }
}

QList<QGraphicsItem*>
Scene::getInteractiveItems(const InteractiveItemFlagType types) const
{
    // We expect all objects that inherit from AbstractGraphicsItem to
    // be interactive (atoms, bonds, etc.), and all objects that don't
    // to be purely graphical (selection highlighting paths, etc.)
    QList<QGraphicsItem*> interactive_items;
    for (auto item : m_interactive_items) {
        if (item_matches_type_flag(item, types)) {
            interactive_items.append(item);
        }
    }
    return interactive_items;
}

void Scene::clearInteractiveItems(const InteractiveItemFlagType types)
{
    // remove all interactive items and reset the rdkit molecule; this will
    // preserve items include selection paths, highlighting items, and
    // potentially the watermark managed by the SketcherWidget
    for (auto* item : getInteractiveItems(types)) {
        removeItem(item);
        m_interactive_items.erase(item);
        delete item;
    }
    if (types & InteractiveItemFlag::ATOM) {
        m_atom_to_atom_item.clear();
    } else if (types & InteractiveItemFlag::BOND) {
        m_bond_to_bond_item.clear();
    } else if (types & InteractiveItemFlag::NON_MOLECULAR) {
        m_non_molecular_to_non_molecular_item.clear();
    }
}

void Scene::updateAtomAndBondItems()
{
    // we need to update all of the atom items before we update any bond
    // items, since bond items pull information from their associated atom
    // items
    for (auto item : getInteractiveItems(InteractiveItemFlag::ATOM)) {
        static_cast<AtomItem*>(item)->updateCachedData();
    }
    updateBondItems();
}

void Scene::updateBondItems()
{
    for (auto item : getInteractiveItems(InteractiveItemFlag::BOND)) {
        static_cast<BondItem*>(item)->updateCachedData();
    }
    // DrawFragmentSceneTool inherits most of the settings from the Scene's
    // AtomItemSettings and BondItemSettings.  Since this method gets called
    // whenever any of those settings are changed, we call updateSceneTool here
    // to regenerate the active scene tool in case DrawFragmentSceneTool is
    // active.  That way, it can inherit the updated settings.
    updateSceneTool();
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
        // check to see whether the mouse has moved far enough to start a
        // drag
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
    if (m_is_during_double_click) {
        // This event is the mouse release for the second click of a double
        // click, so we never got a corresponding mousePressEvent (since it was
        // a mouseDoubleClickEvent instead) and we've already handled the
        // double-click (in mouseDoubleClickEvent).
        m_is_during_double_click = false;
        return;
    }
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
    m_is_during_double_click = true;
    m_left_button_scene_tool->onMouseDoubleClick(event);
}

void Scene::onMouseLeave()
{
    m_left_button_scene_tool->onMouseLeave();
}

void Scene::showContextMenu(QGraphicsSceneMouseEvent* event)
{
    // Collect the set of items that the context menu should interact with
    // based on the position of the cursor
    std::unordered_set<QGraphicsItem*> context_menu_objects;
    auto top_item =
        getTopInteractiveItemAt(event->pos(), InteractiveItemFlag::ALL);
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

AbstractGraphicsItem*
Scene::getTopInteractiveItemAt(const QPointF& pos,
                               const InteractiveItemFlagType types) const
{
    for (auto* item : items(pos)) {
        if (m_interactive_items.count(item) &&
            item_matches_type_flag(item, types)) {
            // all interactive items are AbstractGraphicsItem subclasses, so we
            // can safely use static_cast here
            return static_cast<AbstractGraphicsItem*>(item);
        }
    }
    return nullptr;
}

QList<QGraphicsItem*>
Scene::ensureCompleteAttachmentPoints(const QList<QGraphicsItem*>& items) const
{
    // make sure that .find() calls are O(1) instead of O(N)
    std::unordered_set<QGraphicsItem*> items_set(items.begin(), items.end());
    QList<QGraphicsItem*> new_items;
    for (auto* item : items) {
        if (const auto* atom_item = qgraphicsitem_cast<AtomItem*>(item)) {
            const auto* atom = atom_item->getAtom();
            if (const RDKit::Bond* ap_bond = get_attachment_point_bond(atom)) {
                BondItem* bond_item = m_bond_to_bond_item.at(ap_bond);
                if (items_set.find(bond_item) == items_set.end()) {
                    new_items.push_back(bond_item);
                }
            }
        } else if (auto* bond_item = qgraphicsitem_cast<BondItem*>(item)) {
            const auto* bond = bond_item->getBond();
            if (const RDKit::Atom* ap_atom = get_attachment_point_atom(bond)) {
                AtomItem* atom_item = m_atom_to_atom_item.at(ap_atom);
                if (items_set.find(atom_item) == items_set.end()) {
                    new_items.push_back(atom_item);
                }
            }
        }
    }
    if (new_items.empty()) {
        return items;
    } else {
        return items + new_items;
    }
}

void Scene::onMolModelSelectionChanged()
{
    // block the per-item signals for performance reasons
    Scene::SelectionChangeSignalBlocker signal_blocker(this);
    updateItemSelection();
    updateSelectionHighlighting();
    // signal_blocker emits selectionChanged on destruction
}

void Scene::onModelValuesChanged(const std::unordered_set<ModelKey>& keys)
{
    bool update_atom_and_bond_items = false;

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

            case ModelKey::COLOR_HETEROATOMS:
                setColorHeteroatoms(m_sketcher_model->getValueBool(key));
                update_atom_and_bond_items = true;
                break;
            case ModelKey::SHOW_STEREOCENTER_LABELS:
                m_atom_item_settings.m_stereo_labels_shown =
                    m_sketcher_model->getValueBool(key);
                update_atom_and_bond_items = true;
                break;
            case ModelKey::SHOW_VALENCE_ERRORS:
                m_atom_item_settings.m_valence_errors_shown =
                    m_sketcher_model->getValueBool(key);
                update_atom_and_bond_items = true;
                break;
            default:
                break;
        }
    }
    if (update_atom_and_bond_items) {
        updateAtomAndBondItems();
    }
}

void Scene::setColorHeteroatoms(const bool color_heteroatoms)
{
    if (color_heteroatoms) {
        m_atom_item_settings.setColorScheme(ColorScheme::DEFAULT);
    } else {
        m_atom_item_settings.setMonochromeColorScheme(
            m_atom_item_settings.getAtomColor(-1));
    }
}

void Scene::setFontSize(const qreal size)
{
    m_fonts.setSize(size);
    updateAtomAndBondItems();
}

void Scene::setCarbonsLabeled(const CarbonLabels value)
{
    m_atom_item_settings.m_carbon_labels = value;
    updateAtomAndBondItems();
}

std::shared_ptr<AbstractSceneTool>
Scene::getSceneTool(QGraphicsSceneMouseEvent* const event)
{
    // the button for a release event is stored in event->button(), while
    // the one for a drag event is stored in event->buttons()
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
            case BondTool::SINGLE_OR_DOUBLE:
            case BondTool::SINGLE_OR_AROMATIC:
            case BondTool::DOUBLE_OR_AROMATIC:
            case BondTool::ANY:
                return std::make_shared<DrawBondQuerySceneTool>(bond_tool, this,
                                                                m_mol_model);
            case BondTool::ATOM_CHAIN:
                return std::make_shared<DrawChainSceneTool>(this, m_mol_model);
        }
    } else if (draw_tool == DrawTool::ENUMERATION) {
        switch (m_sketcher_model->getEnumerationTool()) {
            case EnumerationTool::NEW_RGROUP:
                return std::make_shared<DrawIncrementingRGroupSceneTool>(
                    this, m_mol_model);
            case EnumerationTool::EXISTING_RGROUP:
                return std::make_shared<DrawRGroupSceneTool>(
                    m_sketcher_model->getValueInt(ModelKey::RGROUP_NUMBER),
                    this, m_mol_model);
            case EnumerationTool::ATTACHMENT_POINT:
                return std::make_shared<DrawAttachmentPointSceneTool>(
                    this, m_mol_model);
            case EnumerationTool::RXN_ARROW:
                return std::make_shared<ArrowPlusSceneTool>(
                    NonMolecularType::RXN_ARROW, this, m_mol_model);
            case EnumerationTool::RXN_PLUS:
                return std::make_shared<ArrowPlusSceneTool>(
                    NonMolecularType::RXN_PLUS, this, m_mol_model);
            case EnumerationTool::ADD_MAPPING:
                return std::make_shared<AtomMappingSceneTool>(
                    MappingAction::ADD, this, m_mol_model);
            case EnumerationTool::REMOVE_MAPPING:
                return std::make_shared<AtomMappingSceneTool>(
                    MappingAction::REMOVE, this, m_mol_model);
        }
    } else if (draw_tool == DrawTool::CHARGE) {
        auto charge_tool = m_sketcher_model->getChargeTool();
        return std::make_shared<EditChargeSceneTool>(charge_tool, this,
                                                     m_mol_model);
    } else if (draw_tool == DrawTool::MOVE_ROTATE) {
        return std::make_shared<MoveRotateSceneTool>(this, m_mol_model);
    } else if (draw_tool == DrawTool::RING) {
        auto ring_tool = m_sketcher_model->getRingTool();
        return get_draw_fragment_scene_tool(
            ring_tool, m_fonts, m_atom_item_settings, m_bond_item_settings,
            this, m_mol_model);
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

QRectF Scene::getSelectionRect() const
{
    return m_selection_highlighting_item->shape().boundingRect();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/scene.moc"
