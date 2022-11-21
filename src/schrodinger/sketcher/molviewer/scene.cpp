#include "schrodinger/sketcher/molviewer/scene.h"

#include <GraphMol/CoordGen.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <QFont>
#include <QGraphicsSceneMouseEvent>
#include <QMimeData>
#include <QString>
#include <QUrl>

#include "schrodinger/rdkit_extensions/convert.h"

#include "schrodinger/sketcher/error_dialog.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/file_import_export.h"

using schrodinger::rdkit_extensions::Format;

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

Scene::Scene(QObject* parent) : QGraphicsScene(parent)
{
    m_mol = std::make_shared<RDKit::ROMol>(RDKit::ROMol());
}

void Scene::importText(const std::string& text, Format format)
{
    try {
        auto mol = text_to_rdmol(text, format);
        // TODO: deal with chiral flag viz. SHARED-8774
        // TODO: honor existing coordintes if present
        RDKit::CoordGen::addCoords(*mol);
        loadMol(*mol);

    } catch (const std::exception& exc) {
        auto text = QString("Cannot import structure: ") + exc.what();
        show_error_dialog("Import Error", text, window());
    }
}

void Scene::loadMol(const RDKit::ROMol& mol)
{
    auto shared_mol = std::make_shared<RDKit::ROMol>(mol);
    loadMol(shared_mol);
}

void Scene::loadMol(std::shared_ptr<RDKit::ROMol> mol)
{
    clear();
    m_mol = mol;
    const std::size_t num_atoms = m_mol->getNumAtoms();
    std::vector<AtomItem*> atom_items;
    atom_items.reserve(num_atoms);

    // create atom items
    for (std::size_t i = 0; i < num_atoms; ++i) {
        auto const atom = m_mol->getAtomWithIdx(i);
        const auto pos = m_mol->getConformer().getAtomPos(i);
        const auto cur_atom_item =
            new AtomItem(atom, m_fonts, m_atom_item_settings);
        cur_atom_item->setPos(pos.x * VIEW_SCALE, pos.y * VIEW_SCALE);
        addItem(cur_atom_item);
        atom_items.push_back(cur_atom_item);
    }

    // create bond items
    for (auto bond : m_mol->bonds()) {
        const auto* from_atom_item = atom_items[bond->getBeginAtomIdx()];
        const auto* to_atom_item = atom_items[bond->getEndAtomIdx()];
        const auto bond_item = new BondItem(
            bond, *from_atom_item, *to_atom_item, m_bond_item_settings);
        addItem(bond_item);
    }
}

std::shared_ptr<RDKit::ROMol> Scene::getRDKitMolecule() const
{
    return m_mol;
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
            importText(contents, format);
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
    event->accept();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/scene.moc"
