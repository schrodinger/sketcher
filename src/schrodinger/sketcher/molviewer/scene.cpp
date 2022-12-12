#include "schrodinger/sketcher/molviewer/scene.h"

#include <GraphMol/CoordGen.h>
#include <GraphMol/ROMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

#include <QApplication>
#include <QClipboard>
#include <QFont>
#include <QGraphicsSceneMouseEvent>
#include <QMimeData>
#include <QMimeData>
#include <QString>
#include <QUrl>

#include "schrodinger/rdkit_extensions/convert.h"

#include "schrodinger/sketcher/dialog/error_dialog.h"
#include "schrodinger/sketcher/dialog/file_export_dialog.h"
#include "schrodinger/sketcher/dialog/file_save_image_dialog.h"
#include "schrodinger/sketcher/file_import_export.h"
#include "schrodinger/sketcher/image_generation.h"
#include "schrodinger/sketcher/molviewer/atom_item.h"
#include "schrodinger/sketcher/molviewer/bond_item.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/sketcher_model.h"

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
    clear();
}

void Scene::setModel(SketcherModel* model)
{
    if (m_sketcher_model != nullptr) {
        throw std::runtime_error("The model has already been set");
    }
    m_sketcher_model = model;

    connect(this, &Scene::changed, m_sketcher_model,
            &SketcherModel::sceneContentsChanged);
    connect(m_sketcher_model, &SketcherModel::sceneContentsRequested, this,
            [this]() { return items(); });
}

void Scene::loadMol(const RDKit::ROMol& mol)
{
    auto shared_mol = std::make_shared<RDKit::ROMol>(mol);
    loadMol(shared_mol);
}

void Scene::loadMol(std::shared_ptr<RDKit::ROMol> mol)
{
    // TODO: instead of always clearing, handle adding additional mols
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

void Scene::importText(const std::string& text, Format format)
{
    auto show_import_failure = [&](const auto& exc) {
        auto text = QString("Import Failed: ") + exc.what();
        show_error_dialog("Import Error", text, window());
    };

    boost::shared_ptr<RDKit::RWMol> mol{nullptr};
    try {
        mol = rdkit_extensions::text_to_rdmol(text, format);
    } catch (const std::invalid_argument& exc) {
        show_import_failure(exc);
        return;
    } catch (const std::runtime_error& exc) {
        show_import_failure(exc);
        return;
    }

    // TODO: deal with chiral flag viz. SHARED-8774
    // TODO: honor existing coordinates if present
    RDKit::CoordGen::addCoords(*mol);
    loadMol(*mol);
}

std::string Scene::exportText(Format format)
{
    // TODO: handle reactions
    // TODO: handle toggling between selection/everything
    // TODO: handle export of selection as atom/bond properties
    return rdkit_extensions::rdmol_to_text(*m_mol, format);
}

void Scene::clear()
{
    QGraphicsScene::clear();
    m_mol = std::make_shared<RDKit::ROMol>(RDKit::ROMol());
}

void Scene::onImportTextRequested(const std::string& text, Format format)
{
    if (m_sketcher_model->getNewStructuresReplaceContent()) {
        clear();
    }
    importText(text, format);
}

void Scene::showFileExportDialog()
{
    auto dialog = new FileExportDialog(m_sketcher_model, window());
    connect(dialog, &FileExportDialog::exportTextRequested, this,
            [this](Format format) {
                return QString::fromStdString(exportText(format));
            });
    dialog->show();
}

void Scene::showFileSaveImageDialog()
{
    auto dialog = new FileSaveImageDialog(m_sketcher_model, window());
    connect(dialog, &FileSaveImageDialog::exportImageRequested, this,
            [this](auto format, const auto& opts) {
                // TODO: this call uses the existing sketcherScene class;
                // update to render directly from this Scene instance
                return get_image_bytes(*m_mol, format, opts);
            });
    dialog->show();
}

void Scene::onPasteRequested()
{
    auto data = QApplication::clipboard()->mimeData();
    if (data->hasText()) {
        importText(data->text().toStdString(), Format::AUTO_DETECT);
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
    event->acceptProposedAction();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/scene.moc"
