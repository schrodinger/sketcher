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
// Qt's foreach macro conflicts with Boost's foreach, so disable Qt's here
#undef foreach
#include <boost/foreach.hpp>
#include <boost/range/combine.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
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
        // we don't need to update the metadata since undo_mol was stored with
        // valid metadata
        emit moleculeChanged();
        if (selection_changed) {
            emit selectionChanged();
        }
    };
    doCommand(redo, undo, description);
}

void MolModel::finalizeMoleculeChange(bool selection_changed)
{
    updateMoleculeMetadata();
    emit moleculeChanged();
    if (selection_changed) {
        emit selectionChanged();
    }
}

void MolModel::updateMoleculeMetadata()
{
    m_mol.updatePropertyCache(false);
    RDKit::MolOps::fastFindRings(m_mol);
}

std::vector<int> MolModel::getNextNTags(const size_t count,
                                        int& tag_counter) const
{
    std::vector<int> tags;
    tags.reserve(count);
    for (unsigned int i = 0; i < count; i++) {
        tags.push_back(tag_counter++);
    }
    return tags;
}

void MolModel::addAtom(const Element& element, const RDGeom::Point3D& coords,
                       const RDKit::Bond::BondType& bond_type,
                       const RDKit::Bond::BondDir& bond_dir,
                       const RDKit::Atom* const bound_to_atom)
{
    addAtomChain(element, {coords}, bond_type, bond_dir, bound_to_atom);
}

void MolModel::addAtomChain(const Element& element,
                            const std::vector<RDGeom::Point3D>& coords,
                            const RDKit::Bond::BondType& bond_type,
                            const RDKit::Bond::BondDir& bond_dir,
                            const RDKit::Atom* const bound_to_atom)
{
    size_t num_atoms = coords.size();
    size_t num_bonds = num_atoms;
    if (bound_to_atom == nullptr) {
        --num_bonds;
    }
    std::vector<int> atom_tags = getNextNTags(num_atoms, m_next_atom_tag);
    std::vector<int> bond_tags = getNextNTags(num_bonds, m_next_bond_tag);
    int bound_to_atom_tag = -1;
    if (bound_to_atom) {
        bound_to_atom_tag = getTagForAtom(bound_to_atom);
    }
    unsigned int atomic_num = static_cast<unsigned int>(element);
    std::string elem_name = atomic_number_to_name(atomic_num);
    QString desc = QString("Add %1").arg(QString::fromStdString(elem_name));
    auto redo = [this, atom_tags, bond_tags, atomic_num, coords, bond_type,
                 bond_dir, bound_to_atom_tag]() {
        addAtomChainFromCommand(atom_tags, bond_tags, atomic_num, coords,
                                bond_type, bond_dir, bound_to_atom_tag);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::addBond(const RDKit::Atom* const start_atom,
                       const RDKit::Atom* const end_atom,
                       const RDKit::Bond::BondType& bond_type,
                       const RDKit::Bond::BondDir& bond_dir)
{
    int bond_tag = m_next_bond_tag++;
    int start_atom_tag = getTagForAtom(start_atom);
    int end_atom_tag = getTagForAtom(end_atom);

    QString desc = QString("Add bond");
    auto redo = [this, bond_tag, start_atom_tag, end_atom_tag, bond_type,
                 bond_dir]() {
        addBondFromCommand(bond_tag, start_atom_tag, end_atom_tag, bond_type,
                           bond_dir);
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
    if (mol.getNumAtoms() == 0) {
        return;
    }
    std::vector<int> atom_tags =
        getNextNTags(mol.getNumAtoms(), m_next_atom_tag);
    std::vector<int> bond_tags =
        getNextNTags(mol.getNumBonds(), m_next_bond_tag);
    auto redo = [this, mol, atom_tags, bond_tags]() {
        addMolFromCommand(mol, atom_tags, bond_tags);
    };
    doCommandWithMolUndo(redo, description);
}

void MolModel::addMolFromText(const std::string& text, Format format)
{
    boost::shared_ptr<RDKit::RWMol> mol{nullptr};
    mol = rdkit_extensions::to_rdkit(text, format);

    // SHARED-8774: Deal with chiral flag
    // SKETCH-1841: Move this rdkit specific logic into sketcher/rdkit
    // Add 2D coordinates only if the molecule does not already have them
    // present (ie specified via molblock, SMILES extension, etc.)
    auto conformer_2d =
        std::find_if_not(mol->beginConformers(), mol->endConformers(),
                         std::mem_fn(&RDKit::Conformer::is3D));
    if (conformer_2d == mol->endConformers()) {
        RDKit::CoordGen::addCoords(*mol);
    }
    assign_CIP_labels(*mol);
    addMol(*mol);
}

void MolModel::mutateAtom(const RDKit::Atom* const atom, const Element& element)
{
    int atom_tag = getTagForAtom(atom);
    unsigned int atomic_num = static_cast<unsigned int>(element);
    auto redo = [this, atom_tag, atomic_num]() {
        mutateAtomFromCommand(atom_tag, atomic_num);
    };
    doCommandWithMolUndo(redo, "Mutate atom");
}

void MolModel::mutateBond(const RDKit::Bond* const bond,
                          const RDKit::Bond::BondType& bond_type,
                          const RDKit::Bond::BondDir& bond_dir)
{
    int bond_tag = getTagForBond(bond);
    auto redo = [this, bond_tag, bond_type, bond_dir]() {
        mutateBondFromCommand(bond_tag, bond_type, bond_dir);
    };
    doCommandWithMolUndo(redo, "Mutate bond");
}

void MolModel::flipBond(const RDKit::Bond* const bond)
{
    int bond_tag = getTagForBond(bond);
    auto redo = [this, bond_tag]() { flipBondFromCommand(bond_tag); };
    doCommandWithMolUndo(redo, "Flip bond");
}

void MolModel::regenerateCoordinates()
{
    auto cmd = [this]() {
        RDKit::CoordGen::addCoords(m_mol);
        assign_CIP_labels(m_mol);
        emit moleculeChanged();
    };
    doCommandWithMolUndo(cmd, "Clean Up Coordinates");
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
    const QString& desc, std::function<void(RDGeom::Point3D&)> function,
    MergeId merge_id)
{
    auto& conf = m_mol.getConformer();
    auto coords = conf.getPositions();

    // apply function to all coordinates
    std::for_each(coords.begin(), coords.end(), function);

    auto atoms = m_mol.atoms();
    std::vector<int> atom_tags(m_mol.getNumAtoms());
    transform(atoms.begin(), atoms.end(), atom_tags.begin(),
              [this](const RDKit::Atom* atom) { return getTagForAtom(atom); });
    auto current_coords = m_mol.getConformer().getPositions();

    if (merge_id == MergeId::NO_MERGE) {
        auto redo = [this, atom_tags, coords]() {
            this->setCoordinates(atom_tags, coords);
        };
        auto undo = [this, atom_tags, current_coords]() {
            this->setCoordinates(atom_tags, current_coords);
        };
        doCommand(redo, undo, desc);
    } else {
        // issue a mergeable command
        typedef std::tuple<std::vector<int>, std::vector<RDGeom::Point3D>,
                           std::vector<RDGeom::Point3D>>
            MergeData;
        auto redo = [this](MergeData data) {
            auto& atom_tags = std::get<0>(data);
            auto& coords = std::get<1>(data);
            this->setCoordinates(atom_tags, coords);
        };
        auto undo = [this](MergeData data) {
            auto& atom_tags = std::get<0>(data);
            auto& current_coords = std::get<2>(data);
            this->setCoordinates(atom_tags, current_coords);
        };
        auto merge_func = [](MergeData this_data, MergeData other_data) {
            std::get<1>(this_data) = std::get<1>(other_data);
            return this_data;
        };
        MergeData command_data =
            std::make_tuple(atom_tags, coords, current_coords);
        doMergeableCommand<MergeData>(redo, undo, merge_func,
                                      static_cast<int>(merge_id), command_data,
                                      desc);
    }
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

void MolModel::rotateByAngle(float angle)
{
    auto centroid = findCentroid();
    auto angle_radians = angle * M_PI / 180.0; // Convert angle to radians
    auto cosine = std::cos(angle_radians);
    auto sine = std::sin(angle_radians);

    auto rotate = [centroid, sine, cosine](auto& coord) {
        coord -= centroid;
        auto rotated_coord = coord;
        rotated_coord.x = coord.x * cosine - coord.y * sine;
        rotated_coord.y = coord.x * sine + coord.y * cosine;
        coord = rotated_coord + centroid;
    };
    transformCoordinatesWithFunction("Rotate", rotate, MergeId::ROTATE);
}

void MolModel::translateByVector(const RDGeom::Point3D& vector)
{

    auto translate = [vector](auto& coord) { coord += vector; };
    transformCoordinatesWithFunction("Translate", translate,
                                     MergeId::TRANSLATE);
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

void MolModel::addAtomChainFromCommand(
    const std::vector<int>& atom_tags, const std::vector<int>& bond_tags,
    const unsigned int atomic_num, const std::vector<RDGeom::Point3D>& coords,
    const RDKit::Bond::BondType& bond_type,
    const RDKit::Bond::BondDir& bond_dir, const int bound_to_atom_tag)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* prev_atom = nullptr;
    if (bound_to_atom_tag >= 0) {
        prev_atom = m_mol.getUniqueAtomWithBookmark(bound_to_atom_tag);
    }

    int cur_atom_tag;
    RDGeom::Point3D cur_coords;
    auto bond_tags_iter = bond_tags.cbegin();
    BOOST_FOREACH (boost::tie(cur_atom_tag, cur_coords),
                   boost::combine(atom_tags, coords)) {
        auto* atom = new RDKit::Atom(atomic_num);
        unsigned int atom_index = m_mol.addAtom(atom, /* updateLabel = */ false,
                                                /* takeOwnership = */ true);
        setTagForAtom(atom, cur_atom_tag);
        m_mol.getConformer().setAtomPos(atom_index, cur_coords);

        if (prev_atom) {
            const unsigned int cur_bond_tag = *(bond_tags_iter++);
            unsigned int bond_index =
                m_mol.addBond(prev_atom, atom, bond_type) - 1;
            RDKit::Bond* bond = m_mol.getBondWithIdx(bond_index);
            bond->setBondDir(bond_dir);
            setTagForBond(bond, cur_bond_tag);
        }
        prev_atom = atom;
    }
    finalizeMoleculeChange();
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
                                  const RDKit::Bond::BondType& bond_type,
                                  const RDKit::Bond::BondDir& bond_dir)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* start_atom = m_mol.getUniqueAtomWithBookmark(start_atom_tag);
    RDKit::Atom* end_atom = m_mol.getUniqueAtomWithBookmark(end_atom_tag);
    // addBond returns the new number of bonds, *not* the index of the newly
    // added bond, so we have to subtract one to get the index
    unsigned int bond_index =
        m_mol.addBond(start_atom, end_atom, bond_type) - 1;
    RDKit::Bond* bond = m_mol.getBondWithIdx(bond_index);
    bond->setBondDir(bond_dir);
    setTagForBond(bond, bond_tag);
    finalizeMoleculeChange();
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
    finalizeMoleculeChange(selection_changed);
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
    finalizeMoleculeChange();
}

void MolModel::mutateAtomFromCommand(const int atom_tag,
                                     const unsigned int atomic_num)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
    atom->setAtomicNum(atomic_num);
    finalizeMoleculeChange();
}

void MolModel::mutateBondFromCommand(const int bond_tag,
                                     const RDKit::Bond::BondType& bond_type,
                                     const RDKit::Bond::BondDir& bond_dir)
{
    Q_ASSERT(m_in_command);
    RDKit::Bond* bond = m_mol.getUniqueBondWithBookmark(bond_tag);
    bond->setBondType(bond_type);
    bond->setBondDir(bond_dir);
    finalizeMoleculeChange();
}

void MolModel::flipBondFromCommand(const int bond_tag)
{
    Q_ASSERT(m_in_command);
    RDKit::Bond* bond = m_mol.getUniqueBondWithBookmark(bond_tag);
    unsigned int orig_begin = bond->getBeginAtomIdx();
    unsigned int orig_end = bond->getEndAtomIdx();
    bond->setEndAtomIdx(orig_begin);
    bond->setBeginAtomIdx(orig_end);
    finalizeMoleculeChange();
}

void MolModel::setCoordinates(const std::vector<int>& atom_tags,
                              const std::vector<RDGeom::Point3D>& coords)
{
    Q_ASSERT(m_in_command);
    if (atom_tags.size() == coords.size()) {

        for (unsigned int i = 0; i < atom_tags.size(); ++i) {
            RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tags[i]);
            m_mol.getConformer().setAtomPos(atom->getIdx(), coords[i]);
        }
    } else {
        throw std::invalid_argument("setCoordinates: atom_tags "
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
    finalizeMoleculeChange(selection_changed);
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
