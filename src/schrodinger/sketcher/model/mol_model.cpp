#include "schrodinger/sketcher/model/mol_model.h"

#include <algorithm>

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
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
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

namespace
{

/**
 * @return a new atom using the specified atomic number
 */
std::shared_ptr<RDKit::Atom> make_new_atom(const unsigned int atomic_num)
{
    return std::make_shared<RDKit::Atom>(atomic_num);
}

/**
 * @return a new query atom using the specified query
 */
std::shared_ptr<RDKit::Atom> make_new_query_atom(
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query)
{
    auto atom = std::make_shared<RDKit::QueryAtom>();
    atom->setQuery(atom_query->copy());
    return atom;
}

/**
 * @return a new bond using the specified settings
 */
std::shared_ptr<RDKit::Bond>
make_new_bond(const RDKit::Bond::BondType& bond_type,
              const RDKit::Bond::BondDir& bond_dir)
{
    auto bond = std::make_shared<RDKit::Bond>(bond_type);
    bond->setBondDir(bond_dir);
    return bond;
};

std::shared_ptr<RDKit::Bond> make_new_single_bond()
{
    return make_new_bond(RDKit::Bond::BondType::SINGLE,
                         RDKit::Bond::BondDir::NONE);
}

/**
 * @return a new query bond using the specified query
 */
std::shared_ptr<RDKit::Bond> make_new_query_bond(
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    auto bond = std::make_shared<RDKit::QueryBond>();
    bond->setQuery(bond_query->copy());
    return bond;
};

} // namespace

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
    RDKit::RWMol undo_mol = m_mol;
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
                       const RDKit::Atom* const bound_to_atom,
                       const RDKit::Bond::BondType& bond_type,
                       const RDKit::Bond::BondDir& bond_dir)
{
    addAtomChain(element, {coords}, bound_to_atom, bond_type, bond_dir);
}

void MolModel::addAtom(
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
    const RDGeom::Point3D& coords, const RDKit::Atom* const bound_to_atom,
    const RDKit::Bond::BondType& bond_type,
    const RDKit::Bond::BondDir& bond_dir)
{
    addAtomChain(atom_query, {coords}, bound_to_atom, bond_type, bond_dir);
}

void MolModel::addAtom(
    const Element& element, const RDGeom::Point3D& coords,
    const RDKit::Atom* const bound_to_atom,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    addAtomChain(element, {coords}, bound_to_atom, bond_query);
}

void MolModel::addAtom(
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
    const RDGeom::Point3D& coords, const RDKit::Atom* const bound_to_atom,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    addAtomChain(atom_query, {coords}, bound_to_atom, bond_query);
}

void MolModel::addRGroup(const unsigned int r_group_num,
                         const RDGeom::Point3D& coords,
                         const RDKit::Atom* const bound_to_atom)
{
    addRGroupChain({r_group_num}, {coords}, bound_to_atom);
}

void MolModel::addAttachmentPoint(const RDGeom::Point3D& coords,
                                  const RDKit::Atom* const bound_to_atom)
{
    int atom_tag = m_next_atom_tag++;
    int bond_tag = m_next_bond_tag++;
    int bound_to_atom_tag = getTagForAtom(bound_to_atom);
    unsigned int ap_num = get_next_attachment_point_number(&m_mol);

    auto redo = [this, atom_tag, bond_tag, coords, bound_to_atom_tag,
                 ap_num]() {
        auto create_atom = std::bind(make_new_attachment_point, ap_num);
        auto create_bond = make_new_single_bond;
        addAtomChainFromCommand({atom_tag}, {bond_tag}, create_atom, {coords},
                                create_bond, bound_to_atom_tag);
    };
    doCommandWithMolUndo(redo, "Add attachment point");
}

void MolModel::addAtomChain(const Element& element,
                            const std::vector<RDGeom::Point3D>& coords,
                            const RDKit::Atom* const bound_to_atom,
                            const RDKit::Bond::BondType& bond_type,
                            const RDKit::Bond::BondDir& bond_dir)
{
    // structured binding doesn't work with lambda captures until C++20, so we
    // use ties here instead to unpack the return values (same with the ties in
    // the other overloads of this method)
    std::vector<int> atom_tags;
    std::vector<int> bond_tags;
    int bound_to_atom_tag;
    unsigned int atomic_num;
    QString desc;
    std::tie(atom_tags, bond_tags, bound_to_atom_tag) =
        getAtomAndBondTagsForAddingAtomChain(coords, bound_to_atom);
    std::tie(atomic_num, desc) = getAddElementInfo(element);
    auto redo = [this, atom_tags, bond_tags, atomic_num, coords, bond_type,
                 bond_dir, bound_to_atom_tag]() {
        auto create_atom = std::bind(make_new_atom, atomic_num);
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        addAtomChainFromCommand(atom_tags, bond_tags, create_atom, coords,
                                create_bond, bound_to_atom_tag);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::addAtomChain(
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
    const std::vector<RDGeom::Point3D>& coords,
    const RDKit::Atom* const bound_to_atom,
    const RDKit::Bond::BondType& bond_type,
    const RDKit::Bond::BondDir& bond_dir)
{
    std::vector<int> atom_tags;
    std::vector<int> bond_tags;
    int bound_to_atom_tag;
    std::tie(atom_tags, bond_tags, bound_to_atom_tag) =
        getAtomAndBondTagsForAddingAtomChain(coords, bound_to_atom);
    QString desc = QString("Add %1").arg(
        QString::fromStdString(atom_query->getTypeLabel()));
    auto redo = [this, atom_tags, bond_tags, atom_query, coords, bond_type,
                 bond_dir, bound_to_atom_tag]() {
        auto create_atom = std::bind(make_new_query_atom, atom_query);
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        addAtomChainFromCommand(atom_tags, bond_tags, create_atom, coords,
                                create_bond, bound_to_atom_tag);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::addAtomChain(
    const Element& element, const std::vector<RDGeom::Point3D>& coords,
    const RDKit::Atom* const bound_to_atom,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    std::vector<int> atom_tags;
    std::vector<int> bond_tags;
    int bound_to_atom_tag;
    unsigned int atomic_num;
    QString desc;
    std::tie(atom_tags, bond_tags, bound_to_atom_tag) =
        getAtomAndBondTagsForAddingAtomChain(coords, bound_to_atom);
    std::tie(atomic_num, desc) = getAddElementInfo(element);
    auto redo = [this, atom_tags, bond_tags, atomic_num, coords, bond_query,
                 bound_to_atom_tag]() {
        auto create_atom = std::bind(make_new_atom, atomic_num);
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        addAtomChainFromCommand(atom_tags, bond_tags, create_atom, coords,
                                create_bond, bound_to_atom_tag);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::addAtomChain(
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
    const std::vector<RDGeom::Point3D>& coords,
    const RDKit::Atom* const bound_to_atom,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    std::vector<int> atom_tags;
    std::vector<int> bond_tags;
    int bound_to_atom_tag;
    std::tie(atom_tags, bond_tags, bound_to_atom_tag) =
        getAtomAndBondTagsForAddingAtomChain(coords, bound_to_atom);
    QString desc = QString("Add %1").arg(
        QString::fromStdString(atom_query->getTypeLabel()));
    auto redo = [this, atom_tags, bond_tags, atom_query, coords, bond_query,
                 bound_to_atom_tag]() {
        auto create_atom = std::bind(make_new_query_atom, atom_query);
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        addAtomChainFromCommand(atom_tags, bond_tags, create_atom, coords,
                                create_bond, bound_to_atom_tag);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::addRGroupChain(const std::vector<unsigned int> r_group_nums,
                              const std::vector<RDGeom::Point3D>& coords,
                              const RDKit::Atom* const bound_to_atom)
{
    std::vector<int> atom_tags;
    std::vector<int> bond_tags;
    int bound_to_atom_tag;
    std::tie(atom_tags, bond_tags, bound_to_atom_tag) =
        getAtomAndBondTagsForAddingAtomChain(coords, bound_to_atom);
    auto redo = [this, atom_tags, bond_tags, r_group_nums, coords,
                 bound_to_atom_tag]() {
        auto r_group_iter = r_group_nums.cbegin();
        auto create_atom = [&r_group_iter]() {
            return make_new_r_group(*r_group_iter++);
        };
        auto create_bond = make_new_single_bond;
        addAtomChainFromCommand(atom_tags, bond_tags, create_atom, coords,
                                create_bond, bound_to_atom_tag);
    };

    doCommandWithMolUndo(redo, "Add R group");
}

std::tuple<std::vector<int>, std::vector<int>, int>
MolModel::getAtomAndBondTagsForAddingAtomChain(
    const std::vector<RDGeom::Point3D>& coords,
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
    return {atom_tags, bond_tags, bound_to_atom_tag};
}

std::pair<unsigned int, QString>
MolModel::getAddElementInfo(const Element& element)
{
    unsigned int atomic_num = static_cast<unsigned int>(element);
    std::string elem_name = atomic_number_to_name(atomic_num);
    QString desc = QString("Add %1").arg(QString::fromStdString(elem_name));
    return {atomic_num, desc};
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
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        addBondFromCommand(bond_tag, start_atom_tag, end_atom_tag, create_bond);
    };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::addBond(
    const RDKit::Atom* const start_atom, const RDKit::Atom* const end_atom,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    int bond_tag = m_next_bond_tag++;
    int start_atom_tag = getTagForAtom(start_atom);
    int end_atom_tag = getTagForAtom(end_atom);

    QString desc = QString("Add query bond");
    auto redo = [this, bond_tag, start_atom_tag, end_atom_tag, bond_query]() {
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        addBondFromCommand(bond_tag, start_atom_tag, end_atom_tag, create_bond);
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
    auto [expanded_atoms, expanded_bonds] =
        ensureCompleteAttachmentPoints(atoms, bonds);

    std::vector<int> atom_tags;
    for (auto cur_atom : expanded_atoms) {
        atom_tags.push_back(getTagForAtom(cur_atom));
    }

    std::vector<std::tuple<int, int, int>> bond_tags;
    for (auto cur_bond : expanded_bonds) {
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

    // Add 2D coordinates only if the molecule does not already have them
    // present (ie specified via molblock, SMILES extension, etc.)

    auto& conformer_2d = get_2d_conformer(*mol);
    rescale(conformer_2d, *mol, DEFAULT_MOLVIEWER_BOND_LENGTH);

    assign_CIP_labels(*mol);

    // SHARED-8774: Deal with chiral flag
    rdkit_extensions::add_enhanced_stereo_to_chiral_atoms(*mol);

    addMol(*mol);
}

void MolModel::mutateAtom(const RDKit::Atom* const atom, const Element& element)
{
    int atom_tag = getTagForAtom(atom);
    unsigned int atomic_num = static_cast<unsigned int>(element);
    auto redo = [this, atom_tag, atomic_num]() {
        auto create_atom = std::bind(make_new_atom, atomic_num);
        mutateAtomFromCommand(atom_tag, create_atom);
    };
    doCommandWithMolUndo(redo, "Mutate atom");
}

void MolModel::mutateAtom(
    const RDKit::Atom* const atom,
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query)
{
    int atom_tag = getTagForAtom(atom);
    auto redo = [this, atom_tag, atom_query]() {
        auto create_atom = std::bind(make_new_query_atom, atom_query);
        mutateAtomFromCommand(atom_tag, create_atom);
    };
    doCommandWithMolUndo(redo, "Mutate atom");
}

void MolModel::mutateRGroup(const RDKit::Atom* const atom,
                            const unsigned int r_group_num)
{
    int atom_tag = getTagForAtom(atom);
    auto redo = [this, atom_tag, r_group_num]() {
        auto create_atom = std::bind(make_new_r_group, r_group_num);
        mutateAtomFromCommand(atom_tag, create_atom);
    };
    doCommandWithMolUndo(redo, "Mutate atom");
}

void MolModel::mutateBond(const RDKit::Bond* const bond,
                          const RDKit::Bond::BondType& bond_type,
                          const RDKit::Bond::BondDir& bond_dir)
{
    int bond_tag = getTagForBond(bond);
    auto redo = [this, bond_tag, bond_type, bond_dir]() {
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        mutateBondFromCommand(bond_tag, create_bond);
    };
    doCommandWithMolUndo(redo, "Mutate bond");
}

void MolModel::mutateBond(
    const RDKit::Bond* const bond,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    int bond_tag = getTagForBond(bond);
    auto redo = [this, bond_tag, bond_query]() {
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        mutateBondFromCommand(bond_tag, create_bond);
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
    MergeId merge_id, std::vector<const RDKit::Atom*> atoms)
{

    // if no atoms are provided, apply to all atoms
    if (atoms.empty()) {
        auto all_atoms = m_mol.atoms();
        atoms.insert(atoms.end(), all_atoms.begin(), all_atoms.end());
    }
    auto& conf = m_mol.getConformer();
    std::vector<RDGeom::Point3D> coords;
    for (auto atom : atoms) {
        coords.push_back(conf.getAtomPos(atom->getIdx()));
    }
    auto current_coords = coords;

    // apply function to all coordinates
    std::for_each(coords.begin(), coords.end(), function);
    std::vector<int> atom_tags(atoms.size());
    transform(atoms.begin(), atoms.end(), atom_tags.begin(),
              [this](const RDKit::Atom* atom) { return getTagForAtom(atom); });

    // if no merge_id is provided, issue a non-mergeable command
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

void MolModel::rotateByAngle(float angle, const RDGeom::Point3D& pivot_point,
                             const std::vector<const RDKit::Atom*>& atoms)
{

    auto rotate = [pivot_point, angle](auto& coord) {
        coord = rotate_point(coord, pivot_point, angle);
    };
    transformCoordinatesWithFunction("Rotate", rotate, MergeId::ROTATE, atoms);
}

void MolModel::translateByVector(const RDGeom::Point3D& vector,
                                 const std::vector<const RDKit::Atom*>& atoms)
{

    auto translate = [vector](auto& coord) { coord += vector; };
    transformCoordinatesWithFunction("Translate", translate, MergeId::TRANSLATE,
                                     atoms);
}

void MolModel::flipAllVertical()
{
    auto center = find_centroid(m_mol);
    auto flip_y = [center](auto& coord) { coord.y = 2 * center.y - coord.y; };
    transformCoordinatesWithFunction("Flip All Vertical", flip_y);
}

void MolModel::flipAllHorizontal()
{
    auto center = find_centroid(m_mol);
    auto flip_x = [center](auto& coord) { coord.x = 2 * center.x - coord.x; };
    transformCoordinatesWithFunction("flip All Horizontal", flip_x);
}

void MolModel::select(const std::unordered_set<const RDKit::Atom*>& atoms,
                      const std::unordered_set<const RDKit::Bond*>& bonds,
                      const SelectMode select_mode)
{
    auto [expanded_atoms, expanded_bonds] =
        ensureCompleteAttachmentPoints(atoms, bonds);
    std::unordered_set<int> atom_tags;
    for (auto cur_atom : expanded_atoms) {
        atom_tags.insert(getTagForAtom(cur_atom));
    }
    std::unordered_set<int> bond_tags;
    for (auto cur_bond : expanded_bonds) {
        bond_tags.insert(getTagForBond(cur_bond));
    }
    selectTags(atom_tags, bond_tags, select_mode);
}

void MolModel::selectTags(const std::unordered_set<int>& atom_tags,
                          const std::unordered_set<int>& bond_tags,
                          const SelectMode select_mode)
{
    // note that we do not ensure that attachment points are complete here.
    // That happens in select, not selectTags
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

std::pair<std::unordered_set<const RDKit::Atom*>,
          std::unordered_set<const RDKit::Bond*>>
MolModel::ensureCompleteAttachmentPoints(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds)
{
    std::unordered_set<const RDKit::Atom*> new_atoms(atoms);
    std::unordered_set<const RDKit::Bond*> new_bonds(bonds);
    for (const auto* atom : atoms) {
        if (const RDKit::Bond* ap_bond = get_attachment_point_bond(atom)) {
            new_bonds.insert(ap_bond);
        }
    }
    for (const auto* bond : bonds) {
        if (const RDKit::Atom* ap_atom = get_attachment_point_atom(bond)) {
            new_atoms.insert(ap_atom);
        }
    }
    return {new_atoms, new_bonds};
}

void MolModel::addAtomChainFromCommand(
    const std::vector<int>& atom_tags, const std::vector<int>& bond_tags,
    const std::function<std::shared_ptr<RDKit::Atom>()> create_atom,
    const std::vector<RDGeom::Point3D>& coords,
    const std::function<std::shared_ptr<RDKit::Bond>()> create_bond,
    const int bound_to_atom_tag)
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
        auto atom_shared = create_atom();
        auto* atom = atom_shared.get();
        unsigned int atom_index =
            m_mol.addAtom(atom, /* updateLabel = */ false);
        atom = m_mol.getAtomWithIdx(atom_index);
        setTagForAtom(atom, cur_atom_tag);
        m_mol.getConformer().setAtomPos(atom_index, cur_coords);

        if (prev_atom != nullptr) {
            const unsigned int cur_bond_tag = *(bond_tags_iter++);
            auto bond_shared = create_bond();
            auto* bond = bond_shared.get();
            bond->setOwningMol(m_mol);
            bond->setBeginAtom(prev_atom);
            bond->setEndAtom(atom);
            unsigned int bond_index = m_mol.addBond(bond) - 1;
            bond = m_mol.getBondWithIdx(bond_index);
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
    bool selection_changed = m_selected_atom_tags.erase(atom_tag);
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

void MolModel::addBondFromCommand(
    const int bond_tag, const int start_atom_tag, const int end_atom_tag,
    const std::function<std::shared_ptr<RDKit::Bond>()> create_bond)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* start_atom = m_mol.getUniqueAtomWithBookmark(start_atom_tag);
    RDKit::Atom* end_atom = m_mol.getUniqueAtomWithBookmark(end_atom_tag);
    auto bond_shared = create_bond();
    auto* bond = bond_shared.get();
    bond->setOwningMol(m_mol);
    bond->setBeginAtom(start_atom);
    bond->setEndAtom(end_atom);
    // addBond returns the new number of bonds, *not* the index of the newly
    // added bond, so we have to subtract one to get the index
    unsigned int bond_index = m_mol.addBond(bond) - 1;
    bond = m_mol.getBondWithIdx(bond_index);
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
    // we have to determine whether we're deleting an attachment point before we
    // delete any of the bonds, since is_attachment_point() will return false
    // for unbound atoms
    bool attachment_point_deleted =
        std::any_of(atom_tags.begin(), atom_tags.end(), [this](int tag) {
            return is_attachment_point(getAtomFromTag(tag));
        });
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
    if (attachment_point_deleted) {
        renumber_attachment_points(&m_mol);
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

    bool attachment_point_added = false;
    for (int cur_atom_tag : atom_tags) {
        RDKit::Atom* atom = m_mol.getAtomWithIdx(atom_index);
        setTagForAtom(atom, cur_atom_tag);
        attachment_point_added =
            attachment_point_added || is_attachment_point(atom);
        ++atom_index;
    }

    for (int cur_bond_tag : bond_tags) {
        RDKit::Bond* bond = m_mol.getBondWithIdx(bond_index);
        setTagForBond(bond, cur_bond_tag);
        ++bond_index;
    }
    if (attachment_point_added) {
        renumber_attachment_points(&m_mol);
    }
    finalizeMoleculeChange();
}

void MolModel::mutateAtomFromCommand(
    const int atom_tag,
    const std::function<std::shared_ptr<RDKit::Atom>()> create_atom)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
    int atom_index = atom->getIdx();
    auto new_atom = create_atom();
    m_mol.replaceAtom(atom_index, new_atom.get());
    // replaceAtom creates a copy, so we need to fetch the "real" new atom
    auto* mutated_atom = m_mol.getAtomWithIdx(atom_index);
    // The bookmark is automatically updated, but the property is not (unless we
    // passed preserveProbs = true, but that overwrites the R-group property)
    mutated_atom->setProp(TAG_PROPERTY, atom_tag);
    finalizeMoleculeChange();
}

void MolModel::mutateBondFromCommand(
    const int bond_tag,
    const std::function<std::shared_ptr<RDKit::Bond>()> create_bond)
{
    Q_ASSERT(m_in_command);
    RDKit::Bond* bond = m_mol.getUniqueBondWithBookmark(bond_tag);
    int bond_index = bond->getIdx();
    auto new_bond = create_bond();
    m_mol.replaceBond(bond_index, new_bond.get());
    // replaceBond creates a copy, so we need to fetch the "real" new bond
    auto* mutated_bond = m_mol.getBondWithIdx(bond_index);
    // The bookmark is automatically updated, but we have to manually copy the
    // bond tag property
    mutated_bond->setProp(TAG_PROPERTY, bond_tag);
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
    emit coordinatesChanged();
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
