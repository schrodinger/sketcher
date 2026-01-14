#include "schrodinger/sketcher/molviewer/scene.h"

#include <functional>
#include <unordered_set>

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <QApplication>
#include <QFont>
#include <QGraphicsSceneMouseEvent>
#include <QMimeData>
#include <QPainter>
#include <QPainterPath>
#include <QPixmap>
#include <QPointF>
#include <QString>
#include <QTransform>
#include <QUndoStack>
#include <QUrl>
#include <QWidget>

#include "schrodinger/rdkit_extensions/file_format.h"

#include "schrodinger/sketcher/dialog/file_import_export.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/molviewer/abstract_graphics_item.h"
#include "schrodinger/sketcher/molviewer/abstract_monomer_item.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/sgroup_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/molviewer/non_molecular_item.h"
#include "schrodinger/sketcher/molviewer/scene_utils.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/tool/arrow_plus_scene_tool.h"
#include "schrodinger/sketcher/tool/atom_mapping_scene_tool.h"
#include "schrodinger/sketcher/tool/attachment_point_scene_tool.h"
#include "schrodinger/sketcher/tool/select_erase_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_atom_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_bond_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_chain_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_fragment_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_monomer_scene_tool.h"
#include "schrodinger/sketcher/tool/draw_r_group_scene_tool.h"
#include "schrodinger/sketcher/tool/edit_charge_scene_tool.h"
#include "schrodinger/sketcher/tool/move_rotate_scene_tool.h"
#include "schrodinger/sketcher/tool/explicit_h_scene_tool.h"
#include "schrodinger/sketcher/molviewer/halo_highlighting_item.h"
#include "schrodinger/sketcher/molviewer/monomer_utils.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"

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
    m_halo_highlighting_item = new QGraphicsItemGroup();
    m_simplified_stereo_label = addText("");
    m_simplified_stereo_label->setDefaultTextColor(Qt::black);
    m_scene_tool = std::make_shared<NullSceneTool>();

    addItem(m_selection_highlighting_item);
    addItem(m_halo_highlighting_item);

    if (m_mol_model == nullptr || m_sketcher_model == nullptr) {
        throw std::runtime_error("Cannot construct without valid models");
    }
    m_fonts.setSize(m_sketcher_model->getFontSize());
    connect(this, &Scene::changed, m_sketcher_model,
            &SketcherModel::interactiveItemsChanged);
    connect(this, &Scene::selectionChanged, m_sketcher_model,
            &SketcherModel::selectionChanged);

    connect(m_mol_model, &MolModel::modelChanged, this, &Scene::updateItems);

    connect(m_mol_model, &MolModel::coordinatesChanged, this,
            &Scene::moveInteractiveItems);

    connect(m_mol_model, &MolModel::selectionChanged, this,
            &Scene::onMolModelSelectionChanged);
    connect(m_sketcher_model, &SketcherModel::backgroundColorChanged, this,
            &Scene::onBackgroundColorChanged);

    connect(m_sketcher_model, &SketcherModel::valuesChanged, this,
            &Scene::onModelValuesChanged);
    connect(m_sketcher_model, &SketcherModel::displaySettingsChanged, this,
            &Scene::onDisplaySettingsChanged);
    connect(m_sketcher_model, &SketcherModel::interactiveItemsRequested, this,
            [this]() { return getInteractiveItems(); });
    connect(m_sketcher_model, &SketcherModel::selectionRequested, this,
            [this]() { return selectedItems(); });
    connect(m_sketcher_model, &SketcherModel::contextMenuAtomsRequested, this,
            [this]() { return m_context_menu_atoms; });

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

void Scene::moveInteractiveItems()
{
    update_conf_for_mol_graphics_items(
        getInteractiveItems(InteractiveItemFlag::ATOM_OR_MONOMER),
        getInteractiveItems(InteractiveItemFlag::BOND_OR_CONNECTOR),
        getInteractiveItems(InteractiveItemFlag::S_GROUP),
        *m_mol_model->getMol());
    for (auto* non_mol_obj : m_mol_model->getNonMolecularObjects()) {
        auto* non_mol_item =
            m_non_molecular_to_non_molecular_item.at(non_mol_obj);
        non_mol_item->setPos(to_scene_xy(non_mol_obj->getCoords()));
        non_mol_item->updateCachedData();
    }
    updateSelectionHighlighting();
    updateHaloHighlighting();
    m_scene_tool->onStructureUpdated();
}

void Scene::updateItems(const WhatChangedType what_changed)
{
    Scene::SelectionChangeSignalBlocker signal_blocker(this);

    if (what_changed & WhatChanged::MOLECULE) {
        updateMonomerLabelSizeOnModel();
        clearInteractiveItems(InteractiveItemFlag::MOLECULAR_OR_MONOMERIC);
        const auto* mol = m_mol_model->getMol();
        std::vector<QGraphicsItem*> all_items;
        auto atom_display_settings_ptr =
            m_sketcher_model->getAtomDisplaySettingsPtr();
        auto bond_display_settings_ptr =
            m_sketcher_model->getBondDisplaySettingsPtr();

        std::tie(all_items, m_atom_to_atom_item, m_bond_to_bond_item,
                 m_bond_to_secondary_connection_item,
                 m_s_group_to_s_group_item) =
            create_graphics_items_for_mol(mol, m_fonts,
                                          *atom_display_settings_ptr,
                                          *bond_display_settings_ptr);

        for (auto* item : all_items) {
            addItem(item);
            m_interactive_items.insert(item);
        }
        clearHovered();
    }
    if (what_changed & WhatChanged::NON_MOL_OBJS) {
        clearInteractiveItems(InteractiveItemFlag::NON_MOLECULAR);
        QColor color =
            m_sketcher_model->getAtomDisplaySettingsPtr()->getAtomColor(-1);
        for (const auto* model_obj : m_mol_model->getNonMolecularObjects()) {
            auto* item = new NonMolecularItem(model_obj, color);
            auto rd_coords = model_obj->getCoords();
            auto qpos = to_scene_xy(rd_coords);
            item->setPos(qpos);
            m_non_molecular_to_non_molecular_item[model_obj] = item;
            addItem(item);
            m_interactive_items.insert(item);
        }
    }
    updateItemSelection();
    updateSelectionHighlighting();
    updateHaloHighlighting();

    m_scene_tool->onStructureUpdated();
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
    for (auto* bond : m_mol_model->getSelectedSecondaryConnections()) {
        m_bond_to_secondary_connection_item[bond]->setSelected(true);
    }
    for (auto* s_group : m_mol_model->getSelectedSGroups()) {
        m_s_group_to_s_group_item[s_group]->setSelected(true);
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

QRectF
Scene::getInteractiveItemsBoundingRect(const InteractiveItemFlagType types,
                                       bool selection_only) const
{
    QRectF bounding_rect;
    const auto items = getInteractiveItems(types);
    for (QGraphicsItem* item : items) {
        if (selection_only && !item->isSelected()) {
            continue;
        }
        bounding_rect |= item->sceneBoundingRect();
    }
    return bounding_rect;
}

bool Scene::isSimplifiedStereoAnnotationVisible() const
{
    return m_simplified_stereo_label->isVisible();
}

QRectF Scene::getSceneItemsBoundingRect() const
{
    auto items_bounding_rect = getInteractiveItemsBoundingRect();
    if (isSimplifiedStereoAnnotationVisible()) {
        items_bounding_rect |= m_simplified_stereo_label->sceneBoundingRect();
    }
    return items_bounding_rect;
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
    if (types & InteractiveItemFlag::ATOM_OR_MONOMER) {
        m_atom_to_atom_item.clear();
    }
    if (types & InteractiveItemFlag::BOND_OR_CONNECTOR) {
        m_bond_to_bond_item.clear();
        m_bond_to_secondary_connection_item.clear();
    }
    if (types & InteractiveItemFlag::S_GROUP) {
        m_s_group_to_s_group_item.clear();
    }
    if (types & InteractiveItemFlag::NON_MOLECULAR) {
        m_non_molecular_to_non_molecular_item.clear();
    }
}

void Scene::onDisplaySettingsChanged()
{
    auto display_settings_ptr = m_sketcher_model->getAtomDisplaySettingsPtr();
    bool show_simplified_stereo_annotation =
        display_settings_ptr->m_show_simplified_stereo_annotation;
    QString simplified_stereo_annotation;
    if (show_simplified_stereo_annotation) {
        std::string note;
        m_mol_model->getMol()->getPropIfPresent(
            RDKit::common_properties::molNote, note);
        simplified_stereo_annotation = QString::fromStdString(note);
    }
    if (simplified_stereo_annotation.toStdString() !=
        display_settings_ptr->m_simplified_stereo_annotation) {
        // save the actual simplified stereo annotation in the display settings,
        // so that each AtomItem can access it and decide whether to show or
        // hide its labels (atomic labels should be hidden only if the
        // annotation is actually shown and not empty)
        auto new_display_settings(*display_settings_ptr);

        new_display_settings.m_simplified_stereo_annotation =
            simplified_stereo_annotation.toStdString();
        // avoid recursive calls to onDisplaySettingsChanged
        QSignalBlocker signal_blocker(m_sketcher_model);
        m_sketcher_model->setAtomDisplaySettings(new_display_settings);
    }

    m_simplified_stereo_label->setPlainText(simplified_stereo_annotation);
    m_simplified_stereo_label->setVisible(
        !simplified_stereo_annotation.isEmpty());
    m_simplified_stereo_label->setPos(
        getInteractiveItemsBoundingRect().bottomLeft());
    m_fonts.setSize(m_sketcher_model->getFontSize());

    // DrawFragmentSceneTool inherits most of the settings from the Scene's
    // AtomDisplaySettings and BondDisplaySettings.  Since this method gets
    // called whenever any of those settings are changed, we call
    // updateSceneTool here to regenerate the active scene tool in case
    // DrawFragmentSceneTool is active.  That way, it can inherit the updated
    // settings.
    updateSceneTool();

    // refresh the scene's items in case something changed: e.g. when valence
    // errors are hidden on a C the label disappears and all bonds need to be
    // updated
    updateItems(WhatChanged::MOLECULE);
}

void Scene::updateSelectionHighlighting()
{
    m_selection_highlighting_item->highlightItems(selectedItems());
}

void Scene::updateSceneToolAfterSelection()
{
    m_scene_tool->onSelectionChanged();
}

void Scene::updateHaloHighlighting()
{
    for (auto* item : m_halo_highlighting_item->childItems()) {
        m_halo_highlighting_item->removeFromGroup(item);
        delete item;
    }
    for (auto [atoms, bonds, color] : m_mol_model->getHaloHighlighting()) {
        // we want two separate items, with different Z values for atoms and
        // bonds, so we can always have atoms drawn on top of bonds of the
        // corresponding color
        auto atom_Z = static_cast<qreal>(ZOrder::ATOM_HIGHLIGHTING);
        auto bond_Z = static_cast<qreal>(ZOrder::BOND_HIGHLIGHTING);
        QList<QGraphicsItem*> atom_items;
        for (auto atom : atoms) {
            atom_items.append(m_atom_to_atom_item.at(atom));
        }
        QList<QGraphicsItem*> bond_items;
        for (auto bond : bonds) {
            bond_items.append(m_bond_to_bond_item.at(bond));
        }
        for (auto [items, Z] : {std::make_pair(atom_items, atom_Z),
                                std::make_pair(bond_items, bond_Z)}) {
            auto* item = new HaloHighlightingItem(Z);
            item->setPen(color);
            item->setBrush(color);
            item->highlightItems(items);
            m_halo_highlighting_item->addToGroup(item);
        }
    };
}

void Scene::mousePressEvent(QGraphicsSceneMouseEvent* event)
{
    m_scene_tool->mousePressEvent(event);
    event->accept();
    m_scene_tool_from_mouse_press = m_scene_tool;
}

void Scene::mouseMoveEvent(QGraphicsSceneMouseEvent* event)
{
    // TODO: don't pass events to the scene tool if a context menu is
    // visible

    m_scene_tool->mouseMoveEvent(event);
    event->accept();
    updateHovered(event->scenePos());
}

void Scene::updateHovered(const QPointF& pos)
{
    // update the hovered atom or bond
    auto* top_graphics_item =
        getTopInteractiveItemAt(pos, InteractiveItemFlag::ALL);
    if (auto* atom_item =
            dynamic_cast<AbstractAtomOrMonomerItem*>(top_graphics_item)) {
        setHovered(atom_item->getAtom());
    } else if (auto* bond_item = dynamic_cast<AbstractBondOrConnectorItem*>(
                   top_graphics_item)) {
        setHovered(bond_item->getBond());
    } else {
        clearHovered();
    }
}

void Scene::setHovered(
    std::variant<std::monostate, const RDKit::Atom*, const RDKit::Bond*>
        new_hovered)
{
    if (new_hovered == m_hovered) {
        // nothing has changed
        return;
    }

    auto old_hovered = m_hovered;
    m_hovered = new_hovered;

    if (std::holds_alternative<const RDKit::Atom*>(old_hovered) &&
        !std::holds_alternative<const RDKit::Atom*>(new_hovered)) {
        // we'd been hovering over an atom, but now we're not
        emit atomHovered(nullptr);
    } else if (std::holds_alternative<const RDKit::Bond*>(old_hovered) &&
               !std::holds_alternative<const RDKit::Bond*>(new_hovered)) {
        // we'd been hovering over an bond, but now we're not
        emit bondHovered(nullptr);
    }

    if (std::holds_alternative<const RDKit::Atom*>(new_hovered)) {
        // we've just started hovering over this atom
        emit atomHovered(std::get<const RDKit::Atom*>(new_hovered));
    } else if (std::holds_alternative<const RDKit::Bond*>(new_hovered)) {
        // we've just started hovering over this bond
        emit bondHovered(std::get<const RDKit::Bond*>(new_hovered));
    }
}

void Scene::clearHovered()
{
    setHovered(std::monostate());
}

void Scene::mouseReleaseEvent(QGraphicsSceneMouseEvent* event)
{
    m_scene_tool->mouseReleaseEvent(event);
    event->accept();
    m_scene_tool_from_mouse_press = nullptr;
}

void Scene::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* event)
{
    m_scene_tool->mouseDoubleClickEvent(event);
    event->accept();
}

void Scene::onMouseLeave()
{
    m_scene_tool->onMouseLeave();
    clearHovered();
}

void Scene::showContextMenu(QGraphicsSceneMouseEvent* event)
{
    auto pos = event->scenePos();
    auto [atoms, bonds, secondary_connections, sgroups, non_molecular_objects] =
        getModelObjects(SceneSubset::HOVERED, &pos);
    emit showContextMenuRequested(event, atoms, bonds, secondary_connections,
                                  sgroups, non_molecular_objects);
}

std::unordered_set<QGraphicsItem*> Scene::getSelectedInteractiveItems() const
{
    std::unordered_set<QGraphicsItem*> selected_items;
    std::copy_if(m_interactive_items.begin(), m_interactive_items.end(),
                 std::inserter(selected_items, selected_items.begin()),
                 [](QGraphicsItem* item) { return item->isSelected(); });
    return selected_items;
}

std::tuple<std::unordered_set<const RDKit::Atom*>,
           std::unordered_set<const RDKit::Bond*>,
           std::unordered_set<const RDKit::Bond*>,
           std::unordered_set<const RDKit::SubstanceGroup*>,
           std::unordered_set<const NonMolecularObject*>>
Scene::getModelObjects(SceneSubset subset, QPointF* pos) const
{
    if ((subset == SceneSubset::HOVERED ||
         subset == SceneSubset::SELECTED_OR_HOVERED) &&
        pos == nullptr) {
        throw std::runtime_error("Hovered items require a QPointF position");
    }

    std::unordered_set<QGraphicsItem*> items;
    switch (subset) {
        case SceneSubset::ALL:
            items = m_interactive_items;
            break;
        case SceneSubset::SELECTION:
            items = getSelectedInteractiveItems();
            break;
        case SceneSubset::HOVERED:
            if (auto top_item =
                    getTopInteractiveItemAt(*pos, InteractiveItemFlag::ALL)) {
                if (top_item->isSelected()) {
                    items = getSelectedInteractiveItems();
                } else {
                    items = {top_item};
                }
            }
            break;
        case SceneSubset::SELECTED_OR_HOVERED:
            items = getSelectedInteractiveItems();
            if (items.empty()) {
                if (auto top_item = getTopInteractiveItemAt(
                        *pos, InteractiveItemFlag::ALL)) {
                    items = {top_item};
                }
            }
            break;
    }

    return get_model_objects_for_graphics_items(items);
}

AbstractGraphicsItem*
Scene::getTopInteractiveItemAt(const QPointF& pos,
                               const InteractiveItemFlagType types) const
{
    for (auto* item : items(pos)) {
        if (m_interactive_items.count(item) &&
            item_matches_type_flag(item, types)) {
            // all interactive items are AbstractGraphicsItem subclasses, so
            // we can safely use static_cast here
            return static_cast<AbstractGraphicsItem*>(item);
        }
    }
    return nullptr;
}

QList<QGraphicsItem*>
Scene::getCollidingItemsUsingBondMidpoints(QGraphicsItem* item) const
{
    auto items = collidingItems(item);
    auto polygon = item->shape().toFillPolygon();

    // Exclude all bonds whose center point is outside the polygon
    auto is_not_bond_with_center_outside_polygon = [&polygon](auto item) {
        auto bond = dynamic_cast<AbstractBondOrConnectorItem*>(item);
        if (bond == nullptr) {
            return true;
        }
        return polygon.containsPoint(bond->mapToScene(bond->getMidpoint()),
                                     Qt::WindingFill);
    };
    QList<QGraphicsItem*> filtered_items;
    std::copy_if(items.begin(), items.end(), std::back_inserter(filtered_items),
                 is_not_bond_with_center_outside_polygon);
    return filtered_items;
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
                auto* bond_item = m_bond_to_bond_item.at(ap_bond);
                if (items_set.find(bond_item) == items_set.end()) {
                    new_items.push_back(bond_item);
                }
            }
        } else if (auto* bond_item = qgraphicsitem_cast<BondItem*>(item)) {
            const auto* bond = bond_item->getBond();
            if (const RDKit::Atom* ap_atom = get_attachment_point_atom(bond)) {
                auto* atom_item = m_atom_to_atom_item.at(ap_atom);
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

void Scene::onBackgroundColorChanged()
{
    bool has_dark_color_scheme = m_sketcher_model->hasDarkColorScheme();
    m_selection_highlighting_item->setPen(has_dark_color_scheme
                                              ? SELECTION_OUTLINE_COLOR_DARK_BG
                                              : SELECTION_OUTLINE_COLOR);
    m_selection_highlighting_item->setBrush(
        has_dark_color_scheme ? SELECTION_BACKGROUND_COLOR_DARK_BG
                              : SELECTION_BACKGROUND_COLOR);
    m_scene_tool->updateColorsAfterBackgroundColorChange(has_dark_color_scheme);
}

void Scene::onMolModelSelectionChanged()
{
    // block the per-item signals for performance reasons
    Scene::SelectionChangeSignalBlocker signal_blocker(this);
    updateItemSelection();
    updateSelectionHighlighting();
    updateSceneToolAfterSelection();
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
            case ModelKey::MONOMER_TOOL_TYPE:
            case ModelKey::AMINO_ACID_TOOL:
            case ModelKey::NUCLEIC_ACID_TOOL:
                updateSceneTool();
                break;
            default:
                break;
        }
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
            return std::make_shared<DrawElementSceneTool>(element, m_fonts,
                                                          this, m_mol_model);
        } else { // atom_tool == AtomTool::AtomQuery
            auto atom_query = m_sketcher_model->getAtomQuery();
            return std::make_shared<DrawAtomQuerySceneTool>(atom_query, m_fonts,
                                                            this, m_mol_model);
        }
    } else if (draw_tool == DrawTool::BOND) {
        auto bond_tool = m_sketcher_model->getBondTool();
        if (BOND_TOOL_BOND_MAP.count(bond_tool)) {
            // a non-query bond
            return std::make_shared<DrawBondSceneTool>(bond_tool, this,
                                                       m_mol_model);
        } else if (BOND_TOOL_QUERY_MAP.count(bond_tool)) {
            // a query bond
            return std::make_shared<DrawBondQuerySceneTool>(bond_tool, m_fonts,
                                                            this, m_mol_model);
        } else if (bond_tool == BondTool::ATOM_CHAIN) {
            return std::make_shared<DrawChainSceneTool>(this, m_mol_model);
        }
    } else if (draw_tool == DrawTool::ENUMERATION) {
        switch (m_sketcher_model->getEnumerationTool()) {
            case EnumerationTool::NEW_RGROUP:
                return std::make_shared<DrawIncrementingRGroupSceneTool>(
                    m_fonts, this, m_mol_model);
            case EnumerationTool::EXISTING_RGROUP:
                return std::make_shared<DrawRGroupSceneTool>(
                    m_sketcher_model->getValueInt(ModelKey::RGROUP_NUMBER),
                    m_fonts, this, m_mol_model);
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
            ring_tool, m_fonts, *m_sketcher_model->getAtomDisplaySettingsPtr(),
            *m_sketcher_model->getBondDisplaySettingsPtr(), this, m_mol_model);
    } else if (draw_tool == DrawTool::EXPLICIT_H) {
        return std::make_shared<ExplicitHsSceneTool>(this, m_mol_model);
    } else if (draw_tool == DrawTool::MONOMER) {
        auto monomer_tool_type = m_sketcher_model->getMonomerToolType();
        if (monomer_tool_type == MonomerToolType::AMINO_ACID) {
            auto tool = m_sketcher_model->getAminoAcidTool();
            auto res_name = AMINO_ACID_TOOL_TO_RES_NAME.at(tool);
            return std::make_shared<DrawMonomerSceneTool>(
                res_name, rdkit_extensions::ChainType::PEPTIDE, m_fonts, this,
                m_mol_model);
        } else {
            auto tool = m_sketcher_model->getNucleicAcidTool();
            if (NUCLEIC_ACID_TOOL_TO_RES_NAME.contains(tool)) {
                // the tool is for a single monomer
                auto res_name = NUCLEIC_ACID_TOOL_TO_RES_NAME.at(tool);
                // HELM considers DNA to be a type of RNA, so we want an RNA
                // chain type regardless of which nucleic acid we're drawing
                return std::make_shared<DrawMonomerSceneTool>(
                    res_name, rdkit_extensions::ChainType::RNA, m_fonts, this,
                    m_mol_model);
            } else {
                // TODO: the tool is for a full nucleotide
            }
        }
    }
    // tool not yet implemented
    return std::make_shared<NullSceneTool>();
}

void Scene::setSceneTool(std::shared_ptr<AbstractSceneTool> new_scene_tool)
{
    // remove any graphics items associated with the old scene tool
    for (auto* item : m_scene_tool->getGraphicsItems()) {
        removeItem(item);
    }
    disconnect(m_scene_tool.get(), &AbstractSceneTool::contextMenuRequested,
               this, &Scene::showContextMenu);
    disconnect(m_scene_tool.get(), &AbstractSceneTool::newCursorHintRequested,
               this, &Scene::newCursorHintRequested);
    disconnect(m_scene_tool.get(), &AbstractSceneTool::atomDragStarted, this,
               &Scene::onAtomDragStarted);
    disconnect(m_scene_tool.get(), &AbstractSceneTool::atomDragFinished, this,
               &Scene::onAtomDragFinished);
    m_scene_tool = new_scene_tool;
    // add graphics items from the new scene tool
    for (auto* item : m_scene_tool->getGraphicsItems()) {
        addItem(item);
    }
    connect(new_scene_tool.get(), &AbstractSceneTool::contextMenuRequested,
            this, &Scene::showContextMenu);
    connect(new_scene_tool.get(), &AbstractSceneTool::newCursorHintRequested,
            this, &Scene::newCursorHintRequested);
    connect(new_scene_tool.get(), &AbstractSceneTool::atomDragStarted, this,
            &Scene::onAtomDragStarted);
    connect(new_scene_tool.get(), &AbstractSceneTool::atomDragFinished, this,
            &Scene::onAtomDragFinished);
    requestCursorHintUpdate();
    // set the correct colors for the new scene tool, but only
    // if this is not a NullSceneTool. This is required so that we
    // don't hit this when destroying the Scene.
    if (dynamic_cast<NullSceneTool*>(m_scene_tool.get()) == nullptr) {
        m_scene_tool->updateColorsAfterBackgroundColorChange(
            m_sketcher_model->hasDarkColorScheme());
    }
}

void Scene::requestCursorHintUpdate()
{
    auto cursor_hint = m_scene_tool->getDefaultCursorPixmap();
    emit newCursorHintRequested(cursor_hint);
}

bool Scene::isDuringAtomDrag()
{
    return m_currently_dragging_atom;
}

void Scene::onAtomDragStarted()
{
    m_currently_dragging_atom = true;
}

void Scene::onAtomDragFinished(const bool were_atoms_merged)
{
    m_currently_dragging_atom = false;
    if (!were_atoms_merged) {
        emit representationChangingAtomDragFinished();
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
            auto file_path = url.toLocalFile().toStdString();
            auto contents = get_file_text(file_path);
            auto format = rdkit_extensions::get_file_format(file_path);
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

void Scene::updateMonomerLabelSizeOnModel()
{
    if (!m_mol_model->isMonomeric()) {
        return;
    }
    std::unordered_map<int, RDGeom::Point3D> sizes;
    for (auto atom : m_mol_model->getMol()->atoms()) {
        if (!is_atom_monomeric(atom)) {
            continue;
        }
        // create a temporary graphics item to figure out the label size
        auto* item = get_monomer_graphics_item(atom, m_fonts);
        const auto bounding_rect = item->boundingRect();
        sizes[atom->getIdx()] =
            RDGeom::Point3D(bounding_rect.width() / VIEW_SCALE,
                            bounding_rect.height() / VIEW_SCALE, 0);
        delete item;
    }
    m_mol_model->setMonomerSizes(sizes);
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/scene.moc"
