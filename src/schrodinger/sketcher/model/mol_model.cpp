#include "schrodinger/sketcher/model/mol_model.h"

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/CoordGen.h>
#include <GraphMol/MolOps.h>
#include <GraphMol/RWMol.h>
#include <QObject>
#include <QUndoStack>
#include <QtGlobal>
#include <boost/range/combine.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/rdkit/stereochemistry.h"

using schrodinger::rdkit_extensions::Format;

namespace schrodinger
{
namespace sketcher
{

/**
 * The name of the RDKit property used to store atom tags and bond tags
 */
const std::string TAG_PROPERTY = "SKETCHER_TAG";

MolModel::MolModel(QUndoStack* const undo_stack, QObject* parent) :
    AbstractUndoableModel(undo_stack, parent)
{
    auto* conf = new RDKit::Conformer();
    // RDKit takes ownership of this conformer, so we don't have to worry about
    // deleting it
    m_mol.addConformer(conf, true);
}

const RDKit::ROMol* MolModel::getMol() const
{
    // TODO: add API to get reactions
    // TODO: add API to get selected mol
    // TODO: add API to export selection as atom/bond properties
    return &m_mol;
}

std::string MolModel::getMolText(Format format)
{
    return rdkit_extensions::to_string(*getMol(), format);
}

std::unordered_set<const RDKit::Atom*> MolModel::getSelectedAtoms() const
{
    std::unordered_set<const RDKit::Atom*> sel_atoms;
    for (int atom_tag : m_selected_atom_tags) {
        sel_atoms.insert(getAtomFromTag(atom_tag));
    }
    return sel_atoms;
}

std::unordered_set<const RDKit::Bond*> MolModel::getSelectedBonds() const
{
    std::unordered_set<const RDKit::Bond*> sel_bonds;
    for (int bond_tag : m_selected_bond_tags) {
        sel_bonds.insert(getBondFromTag(bond_tag));
    }
    return sel_bonds;
}

void MolModel::doCommandWithMolUndo(const std::function<void()> redo,
                                    const QString& description)
{
    RDKit::ROMol undo_mol = m_mol;
    std::unordered_set<int> sel_atom_tags = m_selected_atom_tags;
    std::unordered_set<int> sel_bond_tags = m_selected_bond_tags;
    auto undo = [this, undo_mol, sel_atom_tags, sel_bond_tags]() {
        m_mol = undo_mol;
        bool selection_changed = sel_atom_tags != m_selected_atom_tags ||
                                 sel_bond_tags != m_selected_bond_tags;
        if (selection_changed) {
            m_selected_atom_tags = sel_atom_tags;
            m_selected_bond_tags = sel_bond_tags;
        }
        emit moleculeChanged();
        if (selection_changed) {
            emit selectionChanged();
        }
    };
    doCommand(redo, undo, description);
}

void MolModel::addAtom(const std::string& element,
                       const RDGeom::Point3D& coords)
{
    // TODO: validate element?
    unsigned int atom_tag = m_next_atom_tag++;

    QString desc = QString("Add %1").arg(QString::fromStdString(element));
    auto redo = [this, atom_tag, element, coords]() {
        addAtomFromCommand(atom_tag, element, coords);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::addBond(const RDKit::Atom* const start_atom,
                       const RDKit::Atom* const end_atom,
                       const RDKit::Bond::BondType& bond_type)
{
    int bond_tag = m_next_bond_tag++;
    int start_atom_tag = getTagForAtom(start_atom);
    int end_atom_tag = getTagForAtom(end_atom);

    QString desc = QString("Add bond");
    auto redo = [this, bond_tag, start_atom_tag, end_atom_tag, bond_type]() {
        addBondFromCommand(bond_tag, start_atom_tag, end_atom_tag, bond_type);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::removeAtomsAndBonds(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds)
{
    if (atoms.empty() && bonds.empty()) {
        return;
    }

    std::vector<int> atom_tags;
    for (auto cur_atom : atoms) {
        atom_tags.push_back(getTagForAtom(cur_atom));
    }

    std::vector<std::tuple<int, int, int>> bond_tags;
    for (auto cur_bond : bonds) {
        int bond_tag = getTagForBond(cur_bond);
        int start_atom_tag = getTagForAtom(cur_bond->getBeginAtom());
        int end_atom_tag = getTagForAtom(cur_bond->getEndAtom());
        bond_tags.push_back({bond_tag, start_atom_tag, end_atom_tag});
    }

    QString desc = QString("Erase");
    auto redo = [this, atom_tags, bond_tags]() {
        removeAtomsAndBondsFromCommand(atom_tags, bond_tags);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::addMol(const RDKit::ROMol& mol, const QString& description)
{
    // generate atom tags for all atoms
    unsigned int num_atoms = mol.getNumAtoms();
    std::vector<int> atom_tags;
    atom_tags.reserve(num_atoms);
    for (unsigned int i = 0; i < num_atoms; i++) {
        atom_tags.push_back(m_next_atom_tag++);
    }

    // generate bond tags for all bonds
    unsigned int num_bonds = mol.getNumBonds();
    std::vector<int> bond_tags;
    bond_tags.reserve(num_bonds);
    for (unsigned int i = 0; i < num_bonds; i++) {
        bond_tags.push_back(m_next_bond_tag++);
    }

    auto redo = [this, mol, atom_tags, bond_tags]() {
        addMolFromCommand(mol, atom_tags, bond_tags);
    };
    doCommandWithMolUndo(redo, description);
}

void MolModel::addMolFromText(const std::string& text, Format format)
{
    boost::shared_ptr<RDKit::RWMol> mol{nullptr};
    mol = rdkit_extensions::to_rdkit(text, format);

    // TODO: deal with chiral flag viz. SHARED-8774
    // TODO: honor existing coordinates if present
    RDKit::CoordGen::addCoords(*mol);
    assign_CIP_labels(*mol);
    addMol(*mol);
}

void MolModel::clear()
{
    if (!(m_mol.getNumAtoms() || m_mol.getNumBonds())) {
        // nothing to clear
        return;
    }
    QString desc = QString("Clear");
    auto redo = [this]() { clearFromCommand(); };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::transformCoordinatesWithFunction(
    const QString& desc, std::function<void(RDGeom::Point3D&)> function)
{
    auto& conf = m_mol.getConformer();
    auto coords = conf.getPositions();

    // apply function to all coordinates
    std::for_each(coords.begin(), coords.end(), function);

    auto atoms = m_mol.atoms();
    std::vector<int> atom_tags(m_mol.getNumAtoms());
    transform(atoms.begin(), atoms.end(), atom_tags.begin(),
              [this](const RDKit::Atom* atom) { return getTagForAtom(atom); });

    auto redo = [this, atom_tags, coords]() {
        this->setCoordinatesFromCommand(atom_tags, coords);
    };
    auto current_coords = m_mol.getConformer().getPositions();
    auto undo = [this, atom_tags, current_coords]() {
        this->setCoordinatesFromCommand(atom_tags, current_coords);
    };
    doCommand(redo, undo, desc);
}

RDGeom::Point3D MolModel::findCentroid() const
{
    auto& conf = m_mol.getConformer();
    auto coords = conf.getPositions();

    // Calculate the centroid by averaging the coordinates
    size_t numAtoms = coords.size();
    RDGeom::Point3D centroid;
    for (const auto& coord : coords) {
        centroid += coord;
    }
    // avoid division by zero
    if (numAtoms > 0) {
        centroid /= static_cast<double>(numAtoms);
    }
    return centroid;
}

void MolModel::flipAllVertical()
{
    auto center = findCentroid();
    auto flip_y = [center](auto& coord) { coord.y = 2 * center.y - coord.y; };
    transformCoordinatesWithFunction("Flip All Vertical", flip_y);
}

void MolModel::flipAllHorizontal()
{
    auto center = findCentroid();
    auto flip_x = [center](auto& coord) { coord.x = 2 * center.x - coord.x; };
    transformCoordinatesWithFunction("flip All Horizontal", flip_x);
}

void MolModel::select(const std::unordered_set<const RDKit::Atom*>& atoms,
                      const std::unordered_set<const RDKit::Bond*>& bonds,
                      const SelectMode select_mode)
{
    std::unordered_set<int> atom_tags;
    for (auto cur_atom : atoms) {
        atom_tags.insert(getTagForAtom(cur_atom));
    }
    std::unordered_set<int> bond_tags;
    for (auto cur_bond : bonds) {
        bond_tags.insert(getTagForBond(cur_bond));
    }
    selectTags(atom_tags, bond_tags, select_mode);
}

void MolModel::selectTags(const std::unordered_set<int>& atom_tags,
                          const std::unordered_set<int>& bond_tags,
                          const SelectMode select_mode)
{
    bool no_tags_specified = atom_tags.empty() && bond_tags.empty();
    if (select_mode == SelectMode::SELECT_ONLY) {
        if (no_tags_specified) {
            clearSelection();
        } else {
            auto undo_macro_raii = createUndoMacro("Select only");
            clearSelection();
            doSelectionCommand(atom_tags, bond_tags, true, "Select");
        }
        return;
    }
    if (no_tags_specified) {
        return;
    }

    auto [selected_atom_tags, deselected_atom_tags] =
        divideBySelected(atom_tags, m_selected_atom_tags);
    auto [selected_bond_tags, deselected_bond_tags] =
        divideBySelected(bond_tags, m_selected_bond_tags);
    if (select_mode == SelectMode::SELECT) {
        // if we passed atom_tags and bond_tags to doSelectionCommand
        // instead of deselected_atom_tags and deselected_bond_tags, then
        // undoing the command could deselect atoms and bonds that should've
        // remained selected.
        doSelectionCommand(deselected_atom_tags, deselected_bond_tags, true,
                           "Select");
    } else if (select_mode == SelectMode::DESELECT) {
        // if we passed atom_tags and bond_tags to doSelectionCommand
        // instead of selected_atom_tags and selected_bond_tags, then
        // undoing the command could select atoms and bonds that should've
        // remained deselected.
        doSelectionCommand(selected_atom_tags, selected_bond_tags, false,
                           "Deselect");
    } else { // select_mode == SelectMode::TOGGLE
        auto undo_macro_raii = createUndoMacro("Toggle selection");
        doSelectionCommand(deselected_atom_tags, deselected_bond_tags, true,
                           "Select");
        doSelectionCommand(selected_atom_tags, selected_bond_tags, false,
                           "Deselect");
    }
}

std::pair<std::unordered_set<int>, std::unordered_set<int>>
MolModel::divideBySelected(const std::unordered_set<int>& tags_to_divide,
                           const std::unordered_set<int>& selected_tags)
{
    std::unordered_set<int> selected;
    std::unordered_set<int> deselected;
    for (int cur_tag : tags_to_divide) {
        if (selected_tags.find(cur_tag) != selected_tags.end()) {
            selected.insert(cur_tag);
        } else {
            deselected.insert(cur_tag);
        }
    }
    return {selected, deselected};
}

void MolModel::doSelectionCommand(
    const std::unordered_set<int>& filtered_atom_tags,
    const std::unordered_set<int>& filtered_bond_tags, const bool to_select,
    const QString& description)
{
    if (filtered_atom_tags.empty() && filtered_bond_tags.empty()) {
        // nothing to select or deselect
        return;
    }

    auto redo = [this, filtered_atom_tags, filtered_bond_tags, to_select]() {
        setSelectedFromCommand(filtered_atom_tags, filtered_bond_tags,
                               to_select);
    };
    auto undo = [this, filtered_atom_tags, filtered_bond_tags, to_select]() {
        setSelectedFromCommand(filtered_atom_tags, filtered_bond_tags,
                               !to_select);
    };
    doCommand(redo, undo, description);
}

void MolModel::clearSelection()
{
    if (m_selected_atom_tags.empty() && m_selected_bond_tags.empty()) {
        // nothing to clear
        return;
    }
    std::unordered_set<int> sel_atom_tags = m_selected_atom_tags;
    std::unordered_set<int> sel_bond_tags = m_selected_bond_tags;

    QString desc = QString("Clear selection");
    auto redo = [this]() { clearSelectionFromCommand(); };
    auto undo = [this, sel_atom_tags, sel_bond_tags]() {
        m_selected_atom_tags = sel_atom_tags;
        m_selected_bond_tags = sel_bond_tags;
    };
    doCommand(redo, undo, desc);
}

void MolModel::selectAll()
{
    auto [atom_tags_to_select, bond_tags_to_select] = getAllUnselectedTags();
    doSelectionCommand(atom_tags_to_select, bond_tags_to_select, true,
                       "Select all");
}

std::pair<std::unordered_set<int>, std::unordered_set<int>>
MolModel::getAllUnselectedTags()
{
    std::unordered_set<int> deselected_atoms;
    for (const auto& [atom_tag, atom] : *m_mol.getAtomBookmarks()) {
        if (m_selected_atom_tags.find(atom_tag) == m_selected_atom_tags.end()) {
            deselected_atoms.insert(atom_tag);
        }
    }
    std::unordered_set<int> deselected_bonds;
    for (const auto& [bond_tag, bond] : *m_mol.getBondBookmarks()) {
        if (m_selected_bond_tags.find(bond_tag) == m_selected_bond_tags.end()) {
            deselected_bonds.insert(bond_tag);
        }
    }
    return {deselected_atoms, deselected_bonds};
}

void MolModel::invertSelection()
{
    auto [atom_tags_to_select, bond_tags_to_select] = getAllUnselectedTags();
    auto undo_macro_raii = createUndoMacro("Invert selection");
    doSelectionCommand(m_selected_atom_tags, m_selected_bond_tags, false,
                       "Deselect");
    doSelectionCommand(atom_tags_to_select, bond_tags_to_select, true,
                       "Select");
}

void MolModel::setTagForAtom(RDKit::Atom* const atom, const int atom_tag)
{
    m_mol.setAtomBookmark(atom, atom_tag);
    atom->setProp(TAG_PROPERTY, atom_tag);
}

int MolModel::getTagForAtom(const RDKit::Atom* const atom)
{
    return atom->getProp<int>(TAG_PROPERTY);
}

void MolModel::setTagForBond(RDKit::Bond* const bond, const int bond_tag)
{
    m_mol.setBondBookmark(bond, bond_tag);
    bond->setProp(TAG_PROPERTY, bond_tag);
}

int MolModel::getTagForBond(const RDKit::Bond* const bond)
{
    return bond->getProp<int>(TAG_PROPERTY);
}

const RDKit::Atom* MolModel::getAtomFromTag(int atom_tag) const
{
    // RDKit is missing const versions of bookmark getters, even though it
    // has const atom getters that take an atom index.  To get around this,
    // we use const_cast.  (See SHARED-9673.)
    return const_cast<RDKit::RWMol*>(&m_mol)->getUniqueAtomWithBookmark(
        atom_tag);
}

const RDKit::Bond* MolModel::getBondFromTag(int bond_tag) const
{
    // RDKit is missing const versions of bookmark getters, even though it
    // has const bond getters that take a bond index.  To get around this,
    // we use const_cast.  (See SHARED-9673.)
    return const_cast<RDKit::RWMol*>(&m_mol)->getUniqueBondWithBookmark(
        bond_tag);
}

void MolModel::addAtomFromCommand(const int atom_tag,
                                  const std::string& element,
                                  const RDGeom::Point3D& coords)
{
    Q_ASSERT(m_in_command);
    // RDKit will take ownership of this atom after we call addAtom, so we
    // don't have to worry about deleting it
    auto* atom = new RDKit::Atom(element);
    // addAtom returns the atom index of the newly added atom (*not* the new
    // number of atoms)
    unsigned int atom_index = m_mol.addAtom(atom, /* updateLabel = */ false,
                                            /* takeOwnership = */ true);
    setTagForAtom(atom, atom_tag);
    m_mol.getConformer().setAtomPos(atom_index, coords);
    // we don't need to update the ring info here since an unbound atom
    // can't be part of a ring
    emit moleculeChanged();
}

bool MolModel::removeAtomFromCommand(const int atom_tag)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
    bool selection_changed = false;
    selection_changed = m_selected_atom_tags.erase(atom_tag);
    // RDKit automatically deletes all bonds involving this atom, so we have
    // to remove those from the selection as well
    for (auto cur_bond : m_mol.atomBonds(atom)) {
        int bond_tag = getTagForBond(cur_bond);
        if (m_selected_bond_tags.erase(bond_tag)) {
            selection_changed = true;
        }
    }
    m_mol.removeAtom(atom);
    return selection_changed;
}

void MolModel::addBondFromCommand(const int bond_tag, const int start_atom_tag,
                                  const int end_atom_tag,
                                  const RDKit::Bond::BondType& bond_type)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* start_atom = m_mol.getUniqueAtomWithBookmark(start_atom_tag);
    RDKit::Atom* end_atom = m_mol.getUniqueAtomWithBookmark(end_atom_tag);
    // addBond returns the new number of bonds, *not* the index of the newly
    // added bond, so we have to subtract one to get the index
    unsigned int bond_index =
        m_mol.addBond(start_atom, end_atom, bond_type) - 1;
    RDKit::Bond* bond = m_mol.getBondWithIdx(bond_index);
    setTagForBond(bond, bond_tag);
    // update the ring info
    RDKit::MolOps::fastFindRings(m_mol);
    emit moleculeChanged();
}

bool MolModel::removeBondFromCommand(const int bond_tag,
                                     const int start_atom_tag,
                                     const int end_atom_tag)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* start_atom = m_mol.getUniqueAtomWithBookmark(start_atom_tag);
    RDKit::Atom* end_atom = m_mol.getUniqueAtomWithBookmark(end_atom_tag);
    bool selection_changed = m_selected_bond_tags.erase(bond_tag);
    m_mol.removeBond(start_atom->getIdx(), end_atom->getIdx());
    return selection_changed;
}

void MolModel::removeAtomsAndBondsFromCommand(
    const std::vector<int>& atom_tags,
    const std::vector<std::tuple<int, int, int>>& bond_tags_with_atoms)
{
    Q_ASSERT(m_in_command);
    bool selection_changed = false;
    // remove the bonds first so that they don't get implicitly deleted when we
    // remove an atom
    for (auto [bond_tag, start_atom_tag, end_atom_tag] : bond_tags_with_atoms) {
        if (removeBondFromCommand(bond_tag, start_atom_tag, end_atom_tag)) {
            selection_changed = true;
        }
    }

    for (int cur_atom_tag : atom_tags) {
        if (removeAtomFromCommand(cur_atom_tag)) {
            selection_changed = true;
        }
    }

    // update the ring info
    RDKit::MolOps::fastFindRings(m_mol);
    emit moleculeChanged();
    if (selection_changed) {
        emit selectionChanged();
    }
}

void MolModel::addMolFromCommand(const RDKit::ROMol& mol,
                                 const std::vector<int>& atom_tags,
                                 const std::vector<int>& bond_tags)
{
    Q_ASSERT(m_in_command);
    // get the starting index for the atoms and bonds to be inserted
    unsigned int atom_index = m_mol.getNumAtoms();
    unsigned int bond_index = m_mol.getNumBonds();
    m_mol.insertMol(mol);
    // update the ring info, which wasn't transferred from mol
    RDKit::MolOps::fastFindRings(m_mol);

    for (int cur_atom_tag : atom_tags) {
        RDKit::Atom* atom = m_mol.getAtomWithIdx(atom_index);
        setTagForAtom(atom, cur_atom_tag);
        ++atom_index;
    }

    for (int cur_bond_tag : bond_tags) {
        RDKit::Bond* bond = m_mol.getBondWithIdx(bond_index);
        setTagForBond(bond, cur_bond_tag);
        ++bond_index;
    }
    emit moleculeChanged();
}

void MolModel::setCoordinatesFromCommand(
    const std::vector<int>& atom_tags,
    const std::vector<RDGeom::Point3D>& coords)
{
    Q_ASSERT(m_in_command);
    if (atom_tags.size() == coords.size()) {

        for (unsigned int i = 0; i < atom_tags.size(); ++i) {
            RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tags[i]);
            m_mol.getConformer().setAtomPos(atom->getIdx(), coords[i]);
        }
    } else {
        throw std::invalid_argument("setCoordinatesFromCommand: atom_tags "
                                    "and coords must have the same size");
    }
    emit moleculeChanged();
}

void MolModel::clearFromCommand()
{
    Q_ASSERT(m_in_command);
    m_mol.clear();
    bool selection_changed =
        !(m_selected_atom_tags.empty() && m_selected_bond_tags.empty());
    m_selected_atom_tags.clear();
    m_selected_bond_tags.clear();
    emit moleculeChanged();
    if (selection_changed) {
        emit selectionChanged();
    }
}

void MolModel::setSelectedFromCommand(const std::unordered_set<int>& atom_tags,
                                      const std::unordered_set<int>& bond_tags,
                                      const bool selected)
{
    Q_ASSERT(m_in_command);
    if (selected) {
        for (int cur_atom_tag : atom_tags) {
            m_selected_atom_tags.insert(cur_atom_tag);
        }
        for (int cur_bond_tag : bond_tags) {
            m_selected_bond_tags.insert(cur_bond_tag);
        }
    } else {
        for (int cur_atom_tag : atom_tags) {
            m_selected_atom_tags.erase(cur_atom_tag);
        }
        for (int cur_bond_tag : bond_tags) {
            m_selected_bond_tags.erase(cur_bond_tag);
        }
    }
    emit selectionChanged();
}

void MolModel::clearSelectionFromCommand()
{
    Q_ASSERT(m_in_command);
    bool selection_changed =
        !(m_selected_atom_tags.empty() && m_selected_bond_tags.empty());
    m_selected_atom_tags.clear();
    m_selected_bond_tags.clear();
    if (selection_changed) {
        emit selectionChanged();
    }
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/model/mol_model.moc"
