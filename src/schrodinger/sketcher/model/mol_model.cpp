#include "schrodinger/sketcher/model/mol_model.h"

#include <algorithm>

#include <fmt/format.h>

#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/Conformer.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <QObject>
#include <QUndoStack>
#include <QtGlobal>
#include <boost/range/combine.hpp>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/rdkit_extensions/sgroup.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/atoms_and_bonds.h"
#include "schrodinger/sketcher/rdkit/fragment.h"
#include "schrodinger/sketcher/rdkit/mol_update.h"
#include "schrodinger/sketcher/rdkit/periodic_table.h"
#include "schrodinger/sketcher/rdkit/rgroup.h"
#include "schrodinger/sketcher/rdkit/subset.h"

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

/**
 * Remove all internal MolModel properties and bookmarks from the given
 * molecule.
 */
void strip_mol_model_tags(RDKit::ROMol& mol)
{
    for (auto* atom : mol.atoms()) {
        atom->clearProp(TAG_PROPERTY);
        atom->clearProp(USER_COLOR);
    }
    for (auto* bond : mol.bonds()) {
        bond->clearProp(TAG_PROPERTY);
        bond->clearProp(USER_COLOR);
    }
    mol.clearAllAtomBookmarks();
    mol.clearAllBondBookmarks();

    // remove the tags from all S-groups
    for (auto& s_group : getSubstanceGroups(mol)) {
        s_group.clearProp(TAG_PROPERTY);
    }
}

} // namespace

MolModelSnapshot::MolModelSnapshot(
    const RDKit::RWMol& mol, const std::vector<NonMolecularObject>& pluses,
    const std::optional<NonMolecularObject>& arrow,
    const std::unordered_set<AtomTag>& selected_atom_tags,
    const std::unordered_set<BondTag>& selected_bond_tags,
    const std::unordered_set<SGroupTag>& selected_s_group_tags,
    const std::unordered_set<NonMolecularTag>& selected_non_molecular_tags,
    const std::vector<HighlightingInfo>& highlighting_info) :
    m_mol(mol),
    m_pluses(pluses),
    m_arrow(arrow),
    m_selected_atom_tags(selected_atom_tags),
    m_selected_bond_tags(selected_bond_tags),
    m_selected_s_group_tags(selected_s_group_tags),
    m_selected_non_molecular_tags(selected_non_molecular_tags),
    m_highlighting_info(highlighting_info)
{
}

bool MolModelSnapshot::isSelectionIdentical(const MolModelSnapshot& other)
{
    return m_selected_atom_tags == other.m_selected_atom_tags &&
           m_selected_bond_tags == other.m_selected_bond_tags &&
           m_selected_s_group_tags == other.m_selected_s_group_tags &&
           m_selected_non_molecular_tags == other.m_selected_non_molecular_tags;
}

MolModel::MolModel(QUndoStack* const undo_stack, QObject* parent) :
    AbstractUndoableModel(undo_stack, parent)
{
    initializeMol();
}

void MolModel::initializeMol()
{
    auto* conf = new RDKit::Conformer();
    // If we don't explicitly set the conformer to 2D, molfiles will be exported
    // with the 3D flag set
    conf->set3D(false);
    // RDKit takes ownership of this conformer, so we don't have to worry about
    // deleting it
    m_mol.addConformer(conf, true);
}

const RDKit::ROMol* MolModel::getMol() const
{
    return &m_mol;
}

boost::shared_ptr<RDKit::ROMol> MolModel::getMolForExport() const
{
    // TODO: add API to export selection as atom/bond properties
    auto mol_copy = boost::make_shared<RDKit::ROMol>(m_mol);
    strip_mol_model_tags(*mol_copy);
    return mol_copy;
}

boost::shared_ptr<RDKit::ROMol> MolModel::getSelectedMolForExport()
{
    // Update the selection to make sure both ends of all selected bonds
    // are included in the selection, because RDKit will fail to create
    // bonds which are missing ends (SKETCH-1232).
    std::unordered_set<const RDKit::Atom*> atoms_to_select;
    for (auto bond : getSelectedBonds()) {
        atoms_to_select.insert(bond->getBeginAtom());
        atoms_to_select.insert(bond->getEndAtom());
    }
    select(atoms_to_select, {}, {}, {}, SelectMode::SELECT);

    // Copy the entire molecule, then remove all atoms and bonds which are not
    // selected, as to preserve any underlying features in the selection; note
    // that RDKit will automatically remove bonds attached to deleted atoms
    RDKit::RWMol mol_copy(m_mol);
    auto [atom_tags, bond_tags, s_group_tags, non_molecular_tags] =
        getAllUnselectedTags();
    mol_copy.beginBatchEdit();
    for (auto atom_tag : atom_tags) {
        mol_copy.removeAtom(getAtomFromTag(atom_tag)->getIdx());
    }
    mol_copy.commitBatchEdit();
    strip_mol_model_tags(mol_copy);
    return boost::make_shared<RDKit::ROMol>(mol_copy);
}

boost::shared_ptr<RDKit::ChemicalReaction>
MolModel::getReactionForExport() const
{
    if (!m_arrow.has_value()) {
        throw std::runtime_error("No reaction arrow found.");
    }
    auto reaction = createReaction(/* strip_tags = */ true);
    if (reaction->getReactants().empty() || reaction->getProducts().empty()) {
        throw std::runtime_error("Incomplete reactions cannot be copied.");
    }
    return reaction;
}

boost::shared_ptr<RDKit::ChemicalReaction>
MolModel::createReaction(const bool strip_tags) const
{
    // get the midpoint of the arrow
    auto arrow_x = m_arrow->getCoords().x;
    auto all_mols =
        RDKit::MolOps::getMolFrags(m_mol, /* sanitizeFrags = */ false);
    auto reaction = boost::make_shared<RDKit::ChemicalReaction>();
    // sort the molecules left-to-right so that they'll be in visual order
    // in the reaction object
    std::vector<std::pair<boost::shared_ptr<RDKit::ROMol>, RDGeom::Point3D>>
        reactants_and_centroids;
    std::vector<std::pair<boost::shared_ptr<RDKit::ROMol>, RDGeom::Point3D>>
        products_and_centroids;
    for (auto cur_mol : all_mols) {
        if (strip_tags) {
            strip_mol_model_tags(*cur_mol);
        }
        auto cur_centroid = find_centroid(*cur_mol);
        if (cur_centroid.x <= arrow_x) {
            reactants_and_centroids.emplace_back(cur_mol, cur_centroid);
        } else {
            products_and_centroids.emplace_back(cur_mol, cur_centroid);
        }
    }
    auto centroid_x_less_than =
        [](std::pair<boost::shared_ptr<RDKit::ROMol>, RDGeom::Point3D> a,
           std::pair<boost::shared_ptr<RDKit::ROMol>, RDGeom::Point3D> b) {
            return a.second.x < b.second.x;
        };
    std::sort(reactants_and_centroids.begin(), reactants_and_centroids.end(),
              centroid_x_less_than);
    std::sort(products_and_centroids.begin(), products_and_centroids.end(),
              centroid_x_less_than);
    for (auto [reactant, centroid] : reactants_and_centroids) {
        reaction->addReactantTemplate(reactant);
    }
    for (auto [product, centroid] : products_and_centroids) {
        reaction->addProductTemplate(product);
    }
    return reaction;
}

bool MolModel::isEmpty() const
{
    return !m_mol.getNumAtoms() && m_pluses.empty() && !m_arrow.has_value();
}

bool MolModel::hasSelection() const
{
    return !(m_selected_atom_tags.empty() && m_selected_bond_tags.empty() &&
             m_selected_s_group_tags.empty() &&
             m_selected_non_molecular_tags.empty());
}

std::unordered_set<const RDKit::Atom*> MolModel::getSelectedAtoms() const
{
    std::unordered_set<const RDKit::Atom*> sel_atoms;
    for (auto atom_tag : m_selected_atom_tags) {
        sel_atoms.insert(getAtomFromTag(atom_tag));
    }
    return sel_atoms;
}

std::unordered_set<const RDKit::Bond*> MolModel::getSelectedBonds() const
{
    std::unordered_set<const RDKit::Bond*> sel_bonds;
    for (auto bond_tag : m_selected_bond_tags) {
        sel_bonds.insert(getBondFromTag(bond_tag));
    }
    return sel_bonds;
}

std::unordered_set<const RDKit::SubstanceGroup*>
MolModel::getSelectedSGroups() const
{
    std::unordered_set<const RDKit::SubstanceGroup*> sel_s_groups;
    if (m_selected_s_group_tags.empty()) {
        return sel_s_groups;
    }
    for (auto& cur_s_group : getSubstanceGroups(m_mol)) {
        auto cur_tag = getTagForSGroup(cur_s_group);
        if (m_selected_s_group_tags.count(cur_tag)) {
            sel_s_groups.insert(&cur_s_group);
        }
    }
    return sel_s_groups;
}

std::unordered_set<const NonMolecularObject*>
MolModel::getSelectedNonMolecularObjects() const
{
    std::unordered_set<const NonMolecularObject*> selected;

    for (auto tag : m_selected_non_molecular_tags) {
        selected.insert(m_tag_to_non_molecular_object.at(tag));
    }
    return selected;
}

std::vector<std::tuple<std::unordered_set<const RDKit::Atom*>,
                       std::unordered_set<const RDKit::Bond*>, QColor>>
MolModel::getHaloHighlighting() const
{
    std::vector<std::tuple<std::unordered_set<const RDKit::Atom*>,
                           std::unordered_set<const RDKit::Bond*>, QColor>>
        highlights;
    for (const auto& info : m_highlighting_info) {
        std::unordered_set<const RDKit::Atom*> atoms;
        std::unordered_set<const RDKit::Bond*> bonds;
        for (auto atom_tag : info.atom_tags) {
            atoms.insert(getAtomFromTag(atom_tag));
        }
        for (auto bond_tag : info.bond_tags) {
            bonds.insert(getBondFromTag(bond_tag));
        }
        highlights.push_back(std::make_tuple(atoms, bonds, info.color));
    }
    return highlights;
}

void MolModel::clearHaloHighlighting()
{
    m_highlighting_info.clear();
}

void MolModel::addHaloHighlighting(const std::unordered_set<int>& atom_indices,
                                   const std::unordered_set<int>& bond_indices,
                                   const QColor& color)
{
    // translate indices to tags
    std::unordered_set<AtomTag> atoms;
    std::unordered_set<BondTag> bonds;
    for (auto atom_idx : atom_indices) {
        atoms.insert(getTagForAtom(m_mol.getAtomWithIdx(atom_idx)));
    }
    for (auto bond_idx : bond_indices) {
        bonds.insert(getTagForBond(m_mol.getBondWithIdx(bond_idx)));
    }

    m_highlighting_info.push_back(HighlightingInfo(atoms, bonds, color));
}

bool MolModel::hasReactionArrow() const
{
    return m_arrow.has_value();
}

const NonMolecularObject* MolModel::getReactionArrow() const
{
    if (m_arrow.has_value()) {
        return &m_arrow.value();
    } else {
        return nullptr;
    }
}

std::unordered_set<const NonMolecularObject*>
MolModel::getNonMolecularObjects() const
{
    std::unordered_set<const NonMolecularObject*> objs;
    for (const auto& plus : m_pluses) {
        objs.insert(&plus);
    }
    if (m_arrow.has_value()) {
        objs.insert(&m_arrow.value());
    }
    return objs;
}

void MolModel::doCommandUsingSnapshots(const std::function<void()> do_func,
                                       const QString& description,
                                       const WhatChangedType to_be_changed)
{
    Q_ASSERT(!m_allow_edits);
    auto undo_snapshot = takeSnapshot();
    m_allow_edits = true;
    do_func();
    if (to_be_changed & WhatChanged::MOLECULE) {
        update_molecule_on_change(m_mol);
    }
    // We don't need to call updateNonMolecularMetadata here since
    // m_tag_to_non_molecular_object (which is what updateNonMolecularMetadata
    // updates) doesn't get saved in the snapshot.
    // (m_tag_to_non_molecular_object uses pointers, which are invalidated by
    // the restoreSnapshot call below, so calling updateNonMolecularMetadata
    // here wouldn't actually do any good.)  updateNonMolecularMetadata is
    // instead called in restoreSnapshot.
    m_allow_edits = false;
    auto redo_snapshot = takeSnapshot();

    bool selection_changed = !undo_snapshot.isSelectionIdentical(redo_snapshot);
    bool arrow_added =
        !undo_snapshot.m_arrow.has_value() && redo_snapshot.m_arrow.has_value();
    bool arrow_removed =
        undo_snapshot.m_arrow.has_value() && !redo_snapshot.m_arrow.has_value();

    auto undo = [this, undo_snapshot, to_be_changed, selection_changed,
                 arrow_removed]() {
        // undoing the removal of an arrow is equivalent to adding an arrow
        restoreSnapshot(undo_snapshot, to_be_changed, selection_changed,
                        arrow_removed);
    };
    auto redo = [this, redo_snapshot, to_be_changed, selection_changed,
                 arrow_added]() {
        restoreSnapshot(redo_snapshot, to_be_changed, selection_changed,
                        arrow_added);
    };
    doCommand(redo, undo, description);
}

MolModelSnapshot MolModel::takeSnapshot() const
{
    return MolModelSnapshot(m_mol, m_pluses, m_arrow, m_selected_atom_tags,
                            m_selected_bond_tags, m_selected_s_group_tags,
                            m_selected_non_molecular_tags, m_highlighting_info);
}

void MolModel::restoreSnapshot(const MolModelSnapshot& snapshot,
                               const WhatChangedType what_changed,
                               const bool selection_changed,
                               const bool arrow_added)
{
    Q_ASSERT(m_allow_edits);
    if (what_changed & WhatChanged::MOLECULE) {
        m_mol = snapshot.m_mol;
        // we don't need to call update_molecule_on_change since all metadata
        // was updated before we took the snapshot
    }
    if (what_changed & WhatChanged::NON_MOL_OBJS) {
        m_pluses = snapshot.m_pluses;
        m_arrow = snapshot.m_arrow;
        updateNonMolecularMetadata();
    }
    if (selection_changed) {
        m_selected_atom_tags = snapshot.m_selected_atom_tags;
        m_selected_bond_tags = snapshot.m_selected_bond_tags;
        m_selected_s_group_tags = snapshot.m_selected_s_group_tags;
        m_selected_non_molecular_tags = snapshot.m_selected_non_molecular_tags;
    }

    if (what_changed & WhatChanged::MOLECULE) {
        emit moleculeChanged();
    }
    if (what_changed & WhatChanged::NON_MOL_OBJS) {
        emit nonMolecularObjectsChanged();
    }
    if (selection_changed) {
        emit selectionChanged();
    }
    if (arrow_added) {
        emit reactionArrowAdded();
    }
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
    auto bound_to_atom_tag = getTagForAtom(bound_to_atom);
    unsigned int ap_num = get_next_attachment_point_number(&m_mol);

    auto cmd_func = [this, coords, bound_to_atom_tag, ap_num]() {
        auto create_atom = std::bind(make_new_attachment_point, ap_num);
        auto create_bond = make_new_single_bond;
        addAtomChainCommandFunc(create_atom, {coords}, create_bond,
                                bound_to_atom_tag);
    };
    doCommandUsingSnapshots(cmd_func, "Add attachment point",
                            WhatChanged::MOLECULE);
}

void MolModel::addAtomChain(const Element& element,
                            const std::vector<RDGeom::Point3D>& coords,
                            const RDKit::Atom* const bound_to_atom,
                            const RDKit::Bond::BondType& bond_type,
                            const RDKit::Bond::BondDir& bond_dir)
{
    auto bound_to_atom_tag =
        getTagForAtom(bound_to_atom, /* allow_null = */ true);
    // structured binding doesn't work with lambda captures until C++20, so we
    // use ties here instead to unpack the return values (same with the ties in
    // the other overloads of this method)
    unsigned int atomic_num;
    QString desc;
    std::tie(atomic_num, desc) = getAddElementInfo(element);
    auto cmd_func = [this, atomic_num, coords, bond_type, bond_dir,
                     bound_to_atom_tag]() {
        auto create_atom = std::bind(make_new_atom, atomic_num);
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        addAtomChainCommandFunc(create_atom, coords, create_bond,
                                bound_to_atom_tag);
    };
    doCommandUsingSnapshots(cmd_func, desc, WhatChanged::MOLECULE);
}

void MolModel::addAtomChain(
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
    const std::vector<RDGeom::Point3D>& coords,
    const RDKit::Atom* const bound_to_atom,
    const RDKit::Bond::BondType& bond_type,
    const RDKit::Bond::BondDir& bond_dir)
{
    auto bound_to_atom_tag =
        getTagForAtom(bound_to_atom, /* allow_null = */ true);
    QString desc = QString("Add %1").arg(
        QString::fromStdString(atom_query->getTypeLabel()));
    auto cmd_func = [this, atom_query, coords, bond_type, bond_dir,
                     bound_to_atom_tag]() {
        auto create_atom = std::bind(make_new_query_atom, atom_query);
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        addAtomChainCommandFunc(create_atom, coords, create_bond,
                                bound_to_atom_tag);
    };
    doCommandUsingSnapshots(cmd_func, desc, WhatChanged::MOLECULE);
}

void MolModel::addAtomChain(
    const Element& element, const std::vector<RDGeom::Point3D>& coords,
    const RDKit::Atom* const bound_to_atom,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    unsigned int atomic_num;
    QString desc;
    auto bound_to_atom_tag =
        getTagForAtom(bound_to_atom, /* allow_null = */ true);
    std::tie(atomic_num, desc) = getAddElementInfo(element);
    auto cmd_func = [this, atomic_num, coords, bond_query,
                     bound_to_atom_tag]() {
        auto create_atom = std::bind(make_new_atom, atomic_num);
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        addAtomChainCommandFunc(create_atom, coords, create_bond,
                                bound_to_atom_tag);
    };
    doCommandUsingSnapshots(cmd_func, desc, WhatChanged::MOLECULE);
}

void MolModel::addAtomChain(
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
    const std::vector<RDGeom::Point3D>& coords,
    const RDKit::Atom* const bound_to_atom,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    auto bound_to_atom_tag =
        getTagForAtom(bound_to_atom, /* allow_null = */ true);
    QString desc = QString("Add %1").arg(
        QString::fromStdString(atom_query->getTypeLabel()));
    auto cmd_func = [this, atom_query, coords, bond_query,
                     bound_to_atom_tag]() {
        auto create_atom = std::bind(make_new_query_atom, atom_query);
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        addAtomChainCommandFunc(create_atom, coords, create_bond,
                                bound_to_atom_tag);
    };
    doCommandUsingSnapshots(cmd_func, desc, WhatChanged::MOLECULE);
}

void MolModel::addRGroupChain(const std::vector<unsigned int> r_group_nums,
                              const std::vector<RDGeom::Point3D>& coords,
                              const RDKit::Atom* const bound_to_atom)
{
    auto bound_to_atom_tag =
        getTagForAtom(bound_to_atom, /* allow_null = */ true);
    auto cmd_func = [this, r_group_nums, coords, bound_to_atom_tag]() {
        auto r_group_iter = r_group_nums.cbegin();
        auto create_atom = [&r_group_iter]() {
            return rdkit_extensions::make_new_r_group(*r_group_iter++);
        };
        auto create_bond = make_new_single_bond;
        addAtomChainCommandFunc(create_atom, coords, create_bond,
                                bound_to_atom_tag);
    };

    doCommandUsingSnapshots(cmd_func, "Add R group", WhatChanged::MOLECULE);
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
    auto start_atom_tag = getTagForAtom(start_atom);
    auto end_atom_tag = getTagForAtom(end_atom);

    QString desc = QString("Add bond");
    auto cmd_func = [this, start_atom_tag, end_atom_tag, bond_type,
                     bond_dir]() {
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        addBondCommandFunc(start_atom_tag, end_atom_tag, create_bond);
    };
    doCommandUsingSnapshots(cmd_func, desc, WhatChanged::MOLECULE);
}

void MolModel::addBond(
    const RDKit::Atom* const start_atom, const RDKit::Atom* const end_atom,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    auto start_atom_tag = getTagForAtom(start_atom);
    auto end_atom_tag = getTagForAtom(end_atom);

    QString desc = QString("Add query bond");
    auto cmd_func = [this, start_atom_tag, end_atom_tag, bond_query]() {
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        addBondCommandFunc(start_atom_tag, end_atom_tag, create_bond);
    };
    doCommandUsingSnapshots(cmd_func, desc, WhatChanged::MOLECULE);
}

void MolModel::addNonMolecularObject(const NonMolecularType& type,
                                     const RDGeom::Point3D& coords)
{
    if (type == NonMolecularType::RXN_ARROW && m_arrow.has_value()) {
        throw std::runtime_error("Only one arrow allowed");
    }
    auto cmd_func = [this, type, coords]() {
        addNonMolecularObjectCommandFunc(type, coords);
    };
    doCommandUsingSnapshots(cmd_func, "Add non-molecular object",
                            WhatChanged::NON_MOL_OBJS);
}

void MolModel::remove(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const std::unordered_set<const RDKit::SubstanceGroup*>& s_groups,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{
    WhatChangedType to_be_changed = WhatChanged::NOTHING;
    if (!(atoms.empty() && bonds.empty() && s_groups.empty())) {
        to_be_changed |= WhatChanged::MOLECULE;
    }
    if (!non_molecular_objects.empty()) {
        to_be_changed |= WhatChanged::NON_MOL_OBJS;
    }
    if (to_be_changed == WhatChanged::NOTHING) {
        // all sets are empty, so there's nothing to remove
        return;
    }
    auto [expanded_atoms, expanded_bonds] =
        ensureCompleteAttachmentPoints(atoms, bonds);

    std::vector<AtomTag> atom_tags;
    for (auto cur_atom : expanded_atoms) {
        atom_tags.push_back(getTagForAtom(cur_atom));
    }

    std::vector<std::tuple<BondTag, AtomTag, AtomTag>> bond_tags;
    for (auto cur_bond : expanded_bonds) {
        auto bond_tag = getTagForBond(cur_bond);
        auto start_atom_tag = getTagForAtom(cur_bond->getBeginAtom());
        auto end_atom_tag = getTagForAtom(cur_bond->getEndAtom());
        bond_tags.push_back({bond_tag, start_atom_tag, end_atom_tag});
    }

    std::vector<NonMolecularTag> non_molecular_tags;
    for (auto* cur_obj : non_molecular_objects) {
        non_molecular_tags.push_back(cur_obj->getTag());
    }

    QString desc = QString("Erase");
    auto cmd_func = [this, atom_tags, bond_tags, s_groups,
                     non_molecular_tags]() {
        removeCommandFunc(atom_tags, bond_tags, s_groups, non_molecular_tags);
    };

    doCommandUsingSnapshots(cmd_func, desc, to_be_changed);
}

void MolModel::addMol(RDKit::RWMol mol, const QString& description,
                      const bool reposition_mol)
{
    if (mol.getNumAtoms() == 0) {
        return;
    }

    // Ensure the newly added mol has coords and necessary stereo information
    prepare_mol(mol);

    if (reposition_mol) {
        if (!m_mol.getNumAtoms()) {
            // if we don't have any existing structure, center the new molecule
            // at the origin
            center_on_origin(mol);
        } else {
            // otherwise, put the new molecule to the right of the existing
            // molecule
            move_molecule_to_the_right_of(mol, m_mol, IMPORT_SPACING);
        }
    }
    auto cmd_func = [this, mol]() { addMolCommandFunc(mol); };
    doCommandUsingSnapshots(cmd_func, description, WhatChanged::MOLECULE);
}

void MolModel::addReaction(const RDKit::ChemicalReaction& reaction)
{
    if (!reaction.getNumReactantTemplates() &&
        !reaction.getNumProductTemplates()) {
        // this reaction is empty
        return;
    }
    if (m_arrow.has_value()) {
        throw std::runtime_error(
            "Sketcher does not support more than one reaction.");
    }
    auto cmd_func = [this, reaction]() { addReactionCommandFunc(reaction); };
    doCommandUsingSnapshots(cmd_func, "Add reaction", WhatChanged::ALL);
}

void MolModel::addFragment(const RDKit::ROMol& fragment_to_add,
                           const RDKit::Atom* const core_start_atom)
{
    QString desc("Add fragment");
    RDKit::RWMol frag(fragment_to_add);
    auto frag_start_atom = prepare_fragment_for_insertion(frag);
    if (core_start_atom == nullptr) {
        // the fragment isn't attached to any existing structure, so we can just
        // add it as is
        addMol(frag, desc, /* reposition_mol = */ false);
        return;
    }

    std::vector<std::pair<unsigned int, AtomFunc>> mutations_to_core_atoms;
    std::vector<std::pair<unsigned int, BondFunc>> mutations_to_core_bonds;
    std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>
        additions_to_core_bonds;
    AtomIdxToFragBondMap core_to_frag_bonds;
    std::tie(mutations_to_core_atoms, mutations_to_core_bonds,
             additions_to_core_bonds, core_to_frag_bonds) =
        get_fragment_addition_info(frag, frag_start_atom, m_mol,
                                   core_start_atom);

    auto cmd_func = [this, frag, mutations_to_core_atoms,
                     mutations_to_core_bonds, additions_to_core_bonds,
                     core_to_frag_bonds]() {
        addFragmentCommandFunc(frag, mutations_to_core_atoms,
                               mutations_to_core_bonds, additions_to_core_bonds,
                               core_to_frag_bonds);
    };
    doCommandUsingSnapshots(cmd_func, desc, WhatChanged::MOLECULE);
}

void MolModel::addFragmentCommandFunc(
    const RDKit::ROMol& fragment,
    const std::vector<std::pair<unsigned int, AtomFunc>>&
        mutations_to_core_atoms,
    const std::vector<std::pair<unsigned int, BondFunc>>&
        mutations_to_core_bonds,
    const std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>&
        additions_to_core_bonds,
    const AtomIdxToFragBondMap& core_to_frag_bonds_by_idx)
{
    addMolCommandFunc(fragment);
    for (const auto& [core_atom_idx, create_atom] : mutations_to_core_atoms) {
        auto core_atom_tag = getTagForAtom(m_mol.getAtomWithIdx(core_atom_idx));
        mutateAtomCommandFunc(core_atom_tag, create_atom);
    }
    for (const auto& [core_bond_idx, create_bond] : mutations_to_core_bonds) {
        auto core_bond_tag = getTagForBond(m_mol.getBondWithIdx(core_bond_idx));
        mutateBondCommandFunc(core_bond_tag, create_bond);
    }
    for (const auto& [begin_atom_idx, end_atom_idx, create_bond] :
         additions_to_core_bonds) {
        auto begin_atom_tag =
            getTagForAtom(m_mol.getAtomWithIdx(begin_atom_idx));
        auto end_atom_tag = getTagForAtom(m_mol.getAtomWithIdx(end_atom_idx));
        addBondCommandFunc(begin_atom_tag, end_atom_tag, create_bond);
    }

    // add bonds between the core (i.e. the previously existing molecule) and
    // the fragment
    for (const auto& [core_atom_idx, frag_atom_idx_to_bond_info] :
         core_to_frag_bonds_by_idx) {
        auto* core_atom = m_mol.getAtomWithIdx(core_atom_idx);
        for (const auto& [frag_atom_idx, bond_info] :
             frag_atom_idx_to_bond_info) {
            const auto& [bond_func, frag_first] = bond_info;
            auto* frag_atom = m_mol.getAtomWithIdx(frag_atom_idx);
            auto begin_atom_tag = getTagForAtom(frag_atom);
            auto end_atom_tag = getTagForAtom(core_atom);
            if (!frag_first) {
                std::swap(begin_atom_tag, end_atom_tag);
            }
            addBondCommandFunc(begin_atom_tag, end_atom_tag, bond_func);
        }
    }
}

void MolModel::setAtomCharge(const RDKit::Atom* const atom, int charge)
{
    auto cmd_func = [this, atom, charge]() {
        getMutableAtom(atom)->setFormalCharge(charge);
    };
    doCommandUsingSnapshots(cmd_func, "Set atom charge", WhatChanged::MOLECULE);
}

void MolModel::setAtomMapping(
    const std::unordered_set<const RDKit::Atom*>& atoms, int mapping)
{
    if (atoms.empty()) {
        return;
    }

    auto redo = [this, atoms, mapping]() {
        Q_ASSERT(m_allow_edits);
        for (auto atom : atoms) {
            getMutableAtom(atom)->setAtomMapNum(mapping);
        }
    };

    doCommandUsingSnapshots(redo, "Set mapping number", WhatChanged::MOLECULE);
}

void MolModel::updateExplicitHs(ExplicitHActions action,
                                std::unordered_set<const RDKit::Atom*> atoms)
{
    bool add_hs;
    if (action == ExplicitHActions::TOGGLE) {
        std::unordered_set<const RDKit::Atom*> atoms_to_check;
        if (!atoms.empty()) {
            atoms_to_check = atoms;
        } else {
            auto all_atoms = m_mol.atoms();
            atoms_to_check.insert(all_atoms.begin(), all_atoms.end());
        }
        add_hs = has_any_implicit_Hs(atoms_to_check);
    } else {
        add_hs = action == ExplicitHActions::ADD;
    }
    if (add_hs) {
        addExplicitHs(atoms);
    } else {
        removeExplicitHs(atoms);
    }
}

void MolModel::addExplicitHs(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    auto cmd_func = [this, atoms]() { addExplicitHsCommandFunc(atoms); };
    doCommandUsingSnapshots(cmd_func, "Add explicit Hs", WhatChanged::MOLECULE);
}

void MolModel::removeExplicitHs(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{

    auto cmd_func = [this, atoms]() { removeExplicitHsCommandFunc(atoms); };
    doCommandUsingSnapshots(cmd_func, "Remove explicit Hs",
                            WhatChanged::MOLECULE);
}

void MolModel::flipBondStereo(std::unordered_set<const RDKit::Bond*> bonds)
{
    std::unordered_set<const RDKit::Bond*> wedges;
    std::unordered_set<const RDKit::Bond*> dashes;
    for (auto& bond : bonds) {
        if (bond->getBondDir() == RDKit::Bond::BEGINWEDGE) {
            wedges.insert(bond);
        }
        if (bond->getBondDir() == RDKit::Bond::BEGINDASH) {
            dashes.insert(bond);
        }
    }
    if (wedges.empty() && dashes.empty()) {
        return;
    }
    auto undo_macro_raii = createUndoMacro("Flip bonds");
    mutateBonds(wedges, BondTool::SINGLE_DOWN);
    mutateBonds(dashes, BondTool::SINGLE_UP);
}

void MolModel::flipBond(const RDKit::Bond* const bond)
{
    auto cmd_func = [this, bond]() { flipBondCommandFunc(bond); };
    doCommandUsingSnapshots(cmd_func, "Flip bond", WhatChanged::MOLECULE);
}

void MolModel::regenerateCoordinates()
{
    auto cmd = [this]() {
        // if this is a reaction, we'll generate new pluses
        m_pluses.clear();
        m_selected_non_molecular_tags.clear();
        if (!hasReactionArrow()) {
            rdkit_extensions::compute2DCoords(m_mol);
        } else {
            regenerateReactionCoordinatesCommandFunc();
        }
        // when coordinates are updated, the wedged/dashed bonds must
        // be updated to match the new coordinates, or else stereochemistry
        // might be inverted if the new coordinates change the parity around
        // any of the chiral centers.
        rdkit_extensions::wedgeMolBonds(m_mol, &m_mol.getConformer());
    };
    doCommandUsingSnapshots(cmd, "Clean Up Coordinates", WhatChanged::ALL);
}

/**
 * @return a list of x coordinates for each atom in the given mol
 */
static std::vector<double>
get_x_coords(const boost::shared_ptr<RDKit::ROMol> mol)
{
    const auto positions = mol->getConformer().getPositions();
    std::vector<double> x_vals;
    std::transform(positions.begin(), positions.end(),
                   std::back_inserter(x_vals), [](auto pos) { return pos.x; });
    return x_vals;
}

/**
 * @return the maximum x coordinate of all atoms in the given mol
 */
static double max_x_of_mol(const boost::shared_ptr<RDKit::ROMol> mol)
{
    auto x_vals = get_x_coords(mol);
    return *std::max_element(x_vals.begin(), x_vals.end());
}

/**
 * @return the minimum x coordinate of all atoms in the given mol
 */
static double min_x_of_mol(const boost::shared_ptr<RDKit::ROMol> mol)
{
    auto x_vals = get_x_coords(mol);
    return *std::min_element(x_vals.begin(), x_vals.end());
}

void MolModel::regenerateReactionCoordinatesCommandFunc()
{
    Q_ASSERT(m_allow_edits);
    auto reaction = createReaction(/* strip_tags = */ false);
    RDDepict::compute2DCoordsForReaction(*reaction,
                                         /* spacing = */ PLUS_SPACING);
    // move all of the products to the right to make room for the arrow (which
    // is wider than a plus)
    const auto PRODUCT_OFFSET = ARROW_SPACING - PLUS_SPACING;
    for (auto product : reaction->getProducts()) {
        auto& conf = product->getConformer();
        for (auto& pos : conf.getPositions()) {
            pos.x += PRODUCT_OFFSET;
        }
    }
    // copy the coordinates from the cleaned up reaction back to our molecule
    auto& mol_model_conf = m_mol.getConformer();
    for (const auto& molecules :
         {reaction->getReactants(), reaction->getProducts()}) {
        for (const auto& rxn_mol : molecules) {
            const auto& rxn_conf = rxn_mol->getConformer();
            for (const auto* rxn_atom : rxn_mol->atoms()) {
                const auto& pos = rxn_conf.getAtomPos(rxn_atom->getIdx());
                auto atom_tag = AtomTag(rxn_atom->getProp<int>(TAG_PROPERTY));
                auto* mol_model_atom = getAtomFromTag(atom_tag);
                mol_model_conf.setAtomPos(mol_model_atom->getIdx(), pos);
            }
        }
    }

    // generate new pluses (the old pluses were already cleared in
    // regenerateCoordinates)
    for (const auto& molecules :
         {reaction->getReactants(), reaction->getProducts()}) {
        if (molecules.empty()) {
            continue;
        }
        // add pluses to the right of everything but the last molecule (since
        // the last reactant gets an arrow to the right of it and the last
        // product doesn't get anything to the right of it)
        for (auto left_mol_it = molecules.begin();
             left_mol_it < (molecules.end() - 1); ++left_mol_it) {
            auto plus_x = max_x_of_mol(*left_mol_it) + PLUS_SPACING / 2.0;
            m_pluses.push_back({NonMolecularType::RXN_PLUS,
                                {plus_x, 0.0, 0.0},
                                ++m_next_non_molecular_tag});
        }
    }

    // re-position the arrow
    double arrow_x = 0.0;
    if (reaction->getNumReactantTemplates() > 0) {
        // move the arrow to the right of the last reactant
        auto last_reactant = reaction->getReactants().back();
        arrow_x = max_x_of_mol(last_reactant) + ARROW_SPACING / 2.0;
    } else if (reaction->getNumProductTemplates() > 0) {
        // there aren't any reactants, so move the arrow to the left of the
        // first product
        auto first_product = reaction->getProducts().front();
        arrow_x = min_x_of_mol(first_product) - ARROW_SPACING / 2.0;
    }
    // if we have neither reactants nor products, then just leave arrow_x at 0.0
    m_arrow->setCoords({arrow_x, 0.0, 0.0});
}

void MolModel::clear()
{
    WhatChangedType to_be_changed = WhatChanged::NOTHING;
    if (m_mol.getNumAtoms()) {
        to_be_changed |= WhatChanged::MOLECULE;
    }
    if (!m_pluses.empty() || m_arrow.has_value()) {
        to_be_changed |= WhatChanged::NON_MOL_OBJS;
    }
    if (to_be_changed == WhatChanged::NOTHING) {
        return;
    }
    QString desc = QString("Clear");
    auto cmd_func = [this]() { clearCommandFunc(); };
    doCommandUsingSnapshots(cmd_func, desc, to_be_changed);
}

void MolModel::transformCoordinatesWithFunction(
    const QString& desc, std::function<void(RDGeom::Point3D&)> function,
    MergeId merge_id, std::unordered_set<const RDKit::Atom*> atoms,
    std::unordered_set<const NonMolecularObject*> non_molecular_objects)
{
    if (isEmpty()) {
        // there aren't any coordinates to transform
        return;
    }

    // if no objects are provided, apply to all objects
    if (atoms.empty() && non_molecular_objects.empty()) {
        auto all_atoms = m_mol.atoms();
        for (auto atom : all_atoms) {
            atoms.insert(atom);
        }
        non_molecular_objects = getNonMolecularObjects();
    }

    // convert to vectors
    std::vector<const RDKit::Atom*> atom_vec(atoms.begin(), atoms.end());
    std::vector<const NonMolecularObject*> non_mol_vec(
        non_molecular_objects.begin(), non_molecular_objects.end());

    // get tags for all objects since we cannot use snapshots for undo().
    std::vector<AtomTag> atom_tags(atom_vec.size());
    transform(atom_vec.begin(), atom_vec.end(), atom_tags.begin(),
              [this](const RDKit::Atom* atom) { return getTagForAtom(atom); });
    std::vector<NonMolecularTag> non_mol_tags(non_mol_vec.size());
    transform(non_mol_vec.begin(), non_mol_vec.end(), non_mol_tags.begin(),
              [](const NonMolecularObject* non_mol_obj) {
                  return non_mol_obj->getTag();
              });

    // get coordinates
    auto& conf = m_mol.getConformer();
    std::vector<RDGeom::Point3D> atom_coords;
    for (auto atom : atom_vec) {
        atom_coords.push_back(conf.getAtomPos(atom->getIdx()));
    }
    auto current_atom_coords = atom_coords;

    std::vector<RDGeom::Point3D> non_mol_coords;
    for (auto* non_mol_obj : non_mol_vec) {
        non_mol_coords.push_back(non_mol_obj->getCoords());
    }
    auto current_non_mol_coords = non_mol_coords;

    // apply function to all coordinates
    std::for_each(atom_coords.begin(), atom_coords.end(), function);
    std::for_each(non_mol_coords.begin(), non_mol_coords.end(), function);

    // if no merge_id is provided, issue a non-mergeable command
    if (merge_id == MergeId::NO_MERGE) {
        auto redo = [this, atom_tags, atom_coords, non_mol_tags,
                     non_mol_coords]() {
            this->setCoordinates(atom_tags, atom_coords, non_mol_tags,
                                 non_mol_coords);
        };
        auto undo = [this, atom_tags, current_atom_coords, non_mol_tags,
                     current_non_mol_coords]() {
            this->setCoordinates(atom_tags, current_atom_coords, non_mol_tags,
                                 current_non_mol_coords);
        };
        doCommand(redo, undo, desc);
    } else {
        // issue a mergeable command
        typedef std::tuple<
            std::vector<AtomTag>, std::vector<RDGeom::Point3D>,
            std::vector<RDGeom::Point3D>, std::vector<NonMolecularTag>,
            std::vector<RDGeom::Point3D>, std::vector<RDGeom::Point3D>>
            MergeData;
        auto redo = [this](MergeData data) {
            auto& [atom_tags, atom_coords, current_atom_coords, non_mol_tags,
                   non_mol_coords, current_non_mol_coords] = data;
            this->setCoordinates(atom_tags, atom_coords, non_mol_tags,
                                 non_mol_coords);
        };
        auto undo = [this](MergeData data) {
            auto& [atom_tags, atom_coords, current_atom_coords, non_mol_tags,
                   non_mol_coords, current_non_mol_coords] = data;
            this->setCoordinates(atom_tags, current_atom_coords, non_mol_tags,
                                 current_non_mol_coords);
        };
        auto merge_func = [](MergeData this_data, MergeData other_data) {
            // update atom_coords
            std::get<1>(this_data) = std::get<1>(other_data);
            // update non_mol_coords
            std::get<4>(this_data) = std::get<4>(other_data);
            return this_data;
        };
        MergeData command_data = std::make_tuple(
            atom_tags, atom_coords, current_atom_coords, non_mol_tags,
            non_mol_coords, current_non_mol_coords);
        doMergeableCommand<MergeData>(redo, undo, merge_func,
                                      static_cast<int>(merge_id), command_data,
                                      desc);
    }
}

void MolModel::rotateByAngle(
    float angle, const RDGeom::Point3D& pivot_point,
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{

    auto rotate = [pivot_point, angle](auto& coord) {
        coord = rotate_point(coord, pivot_point, angle);
    };
    transformCoordinatesWithFunction("Rotate", rotate, MergeId::ROTATE, atoms,
                                     non_molecular_objects);
}

void MolModel::translateByVector(
    const RDGeom::Point3D& vector,
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{

    auto translate = [vector](auto& coord) { coord += vector; };
    transformCoordinatesWithFunction("Translate", translate, MergeId::TRANSLATE,
                                     atoms, non_molecular_objects);
}

void MolModel::flipSubstituent(const RDKit::Bond* const bond)
{
    auto atoms = get_smaller_substituent_atoms(m_mol, *bond);
    if (atoms.empty()) {
        return;
    }
    flipAroundSegment(
        m_mol.getConformer().getAtomPos(bond->getBeginAtom()->getIdx()),
        m_mol.getConformer().getAtomPos(bond->getEndAtom()->getIdx()), atoms);
}

void MolModel::flipSelection()
{
    auto selected_atoms = getSelectedAtoms();
    if (selected_atoms.empty()) {
        return;
    }
    auto mol = getMol();
    // Find all bonds with one selected atom and one unselected atom
    std::unordered_set<const RDKit::Bond*> crossing_bonds;
    for (auto bond : mol->bonds()) {
        if (selected_atoms.count(bond->getBeginAtom()) !=
            selected_atoms.count(bond->getEndAtom())) {
            crossing_bonds.insert(bond);
        }
    }
    // If there is only one, use it as the point
    if (crossing_bonds.size() == 1) {
        auto bond = *crossing_bonds.begin();
        auto start_coord =
            mol->getConformer().getAtomPos(bond->getBeginAtom()->getIdx());
        auto end_coord =
            mol->getConformer().getAtomPos(bond->getEndAtom()->getIdx());
        flipAroundSegment(start_coord, end_coord, selected_atoms);
    } else {
        // throw an exception
        throw std::runtime_error("Exactly one bond must be selected");
    }
}

void MolModel::flipAroundSegment(
    const RDGeom::Point3D& p1, const RDGeom::Point3D& p2,
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    if (atoms.empty()) {
        return;
    }
    auto flip = [p1, p2](auto& coord) { coord = flip_point(coord, p1, p2); };
    auto undo_macro_raii = createUndoMacro("Flip selection");
    transformCoordinatesWithFunction("", flip, MergeId::NO_MERGE, atoms);
    // get all bonds between atoms and flip wedges to dashes and vice-versa
    std::unordered_set<const RDKit::Bond*> bonds;
    auto& mol = (*atoms.begin())->getOwningMol();
    for (auto& atom : atoms) {
        for (auto& bond : mol.atomBonds(atom)) {
            if (atoms.count(bond->getBeginAtom()) &&
                atoms.count(bond->getEndAtom())) {
                bonds.insert(bond);
            }
        }
    }
    flipBondStereo(bonds);
}

void MolModel::flipAllHorizontal()
{
    auto undo_macro_raii = createUndoMacro("Flip all horizontal");
    auto center = find_centroid(m_mol, getNonMolecularObjects());
    auto flip_x = [center](auto& coord) { coord.x = 2 * center.x - coord.x; };
    transformCoordinatesWithFunction("", flip_x);
    flipBondStereo({m_mol.beginBonds(), m_mol.endBonds()});
}

void MolModel::flipAllVertical()
{
    auto undo_macro_raii = createUndoMacro("Flip all vertical");
    auto center = find_centroid(m_mol, getNonMolecularObjects());
    auto flip_y = [center](auto& coord) { coord.y = 2 * center.y - coord.y; };
    transformCoordinatesWithFunction("", flip_y);
    flipBondStereo({m_mol.beginBonds(), m_mol.endBonds()});
}

void MolModel::mergeAtoms(
    const std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>>&
        atoms_to_merge)
{
    auto cmd = [this, atoms_to_merge]() {
        mergeAtomsCommandFunc(atoms_to_merge);
    };
    doCommandUsingSnapshots(cmd, "Merge dragged atoms", WhatChanged::MOLECULE);
}

void MolModel::mergeAtomsCommandFunc(
    const std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>>&
        atoms_to_merge)
{
    auto& conf = m_mol.getConformer();
    for (auto [atom_to_keep, stationary_atom] : atoms_to_merge) {
        auto atom_to_keep_idx = atom_to_keep->getIdx();
        auto stationary_atom_idx = stationary_atom->getIdx();

        // move the atom_to_keep to stationary_atom's coordinates
        auto& coords = conf.getAtomPos(stationary_atom_idx);
        conf.setAtomPos(atom_to_keep_idx, coords);

        // move all bonds of stationary_atom to atom_to_keep
        for (auto* bond : m_mol.atomBonds(stationary_atom)) {
            auto other_atom_idx = bond->getOtherAtomIdx(stationary_atom_idx);
            if (m_mol.getBondBetweenAtoms(atom_to_keep_idx, other_atom_idx)) {
                // There's already a bond between this atom and atom_to_keep, so
                // skip this bond.  It'll get deleted when we remove
                // stationary_atom.
                continue;
            }
            // the molecule won't update its connectivity graph if we edit the
            // bond in place, so we instead remove it and add a copy
            auto* replacement_bond = new RDKit::Bond(*bond);
            if (bond->getBeginAtomIdx() == stationary_atom_idx) {
                replacement_bond->setBeginAtomIdx(atom_to_keep_idx);
            } else {
                replacement_bond->setEndAtomIdx(atom_to_keep_idx);
            }
            auto bond_tag = getTagForBond(bond);
            m_mol.removeBond(bond->getBeginAtomIdx(), bond->getEndAtomIdx());
            m_mol.addBond(replacement_bond, /* take_ownership = */ true);
            setTagForBond(replacement_bond, bond_tag);
        }

        // erase stationary_atom
        removeAtomCommandFunc(getTagForAtom(stationary_atom));
    }
}

void MolModel::select(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const std::unordered_set<const RDKit::SubstanceGroup*>& s_groups,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects,
    const SelectMode select_mode)
{
    auto [expanded_atoms, expanded_bonds] =
        ensureCompleteAttachmentPoints(atoms, bonds);
    std::unordered_set<AtomTag> atom_tags;
    for (auto* cur_atom : expanded_atoms) {
        atom_tags.insert(getTagForAtom(cur_atom));
    }
    std::unordered_set<BondTag> bond_tags;
    for (auto* cur_bond : expanded_bonds) {
        bond_tags.insert(getTagForBond(cur_bond));
    }
    std::unordered_set<SGroupTag> s_group_tags;
    for (auto* cur_s_group : s_groups) {
        s_group_tags.insert(getTagForSGroup(*cur_s_group));
    }
    std::unordered_set<NonMolecularTag> non_molecular_tags;
    for (auto* cur_non_molecular_obj : non_molecular_objects) {
        non_molecular_tags.insert(cur_non_molecular_obj->getTag());
    }
    selectTags(atom_tags, bond_tags, s_group_tags, non_molecular_tags,
               select_mode);
}

void MolModel::selectTags(
    const std::unordered_set<AtomTag>& atom_tags,
    const std::unordered_set<BondTag>& bond_tags,
    const std::unordered_set<SGroupTag>& s_group_tags,
    const std::unordered_set<NonMolecularTag>& non_molecular_tags,
    const SelectMode select_mode)
{
    // note that we do not ensure that attachment points are complete here.
    // That happens in select, not selectTags
    bool no_tags_specified = atom_tags.empty() && bond_tags.empty() &&
                             s_group_tags.empty() && non_molecular_tags.empty();
    if (select_mode == SelectMode::SELECT_ONLY) {
        if (no_tags_specified) {
            clearSelection();
        } else {
            auto undo_macro_raii = createUndoMacro("Select only");
            clearSelection();
            doSelectionCommand(atom_tags, bond_tags, s_group_tags,
                               non_molecular_tags, true, "Select");
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
    auto [selected_s_group_tags, deselected_s_group_tags] =
        divideBySelected(s_group_tags, m_selected_s_group_tags);
    auto [selected_non_molecular_tags, deselected_non_molecular_tags] =
        divideBySelected(non_molecular_tags, m_selected_non_molecular_tags);
    if (select_mode == SelectMode::SELECT) {
        // if we passed atom_tags and bond_tags to doSelectionCommand
        // instead of deselected_atom_tags and deselected_bond_tags, then
        // undoing the command could deselect atoms and bonds that should've
        // remained selected.
        doSelectionCommand(deselected_atom_tags, deselected_bond_tags,
                           deselected_s_group_tags,
                           deselected_non_molecular_tags, true, "Select");
    } else if (select_mode == SelectMode::DESELECT) {
        // if we passed atom_tags and bond_tags to doSelectionCommand
        // instead of selected_atom_tags and selected_bond_tags, then
        // undoing the command could select atoms and bonds that should've
        // remained deselected.
        doSelectionCommand(selected_atom_tags, selected_bond_tags,
                           selected_s_group_tags, selected_non_molecular_tags,
                           false, "Deselect");
    } else { // select_mode == SelectMode::TOGGLE
        auto undo_macro_raii = createUndoMacro("Toggle selection");
        doSelectionCommand(deselected_atom_tags, deselected_bond_tags,
                           deselected_s_group_tags,
                           deselected_non_molecular_tags, true, "Select");
        doSelectionCommand(selected_atom_tags, selected_bond_tags,
                           selected_s_group_tags, selected_non_molecular_tags,
                           false, "Deselect");
    }
}

template <class T> std::pair<std::unordered_set<T>, std::unordered_set<T>>
MolModel::divideBySelected(const std::unordered_set<T>& tags_to_divide,
                           const std::unordered_set<T>& selected_tags)
{
    std::unordered_set<T> selected;
    std::unordered_set<T> deselected;
    for (T cur_tag : tags_to_divide) {
        if (selected_tags.find(cur_tag) != selected_tags.end()) {
            selected.insert(cur_tag);
        } else {
            deselected.insert(cur_tag);
        }
    }
    return {selected, deselected};
}

void MolModel::doSelectionCommand(
    const std::unordered_set<AtomTag>& filtered_atom_tags,
    const std::unordered_set<BondTag>& filtered_bond_tags,
    const std::unordered_set<SGroupTag>& filtered_s_group_tags,
    const std::unordered_set<NonMolecularTag>& filtered_non_molecular_tags,
    const bool to_select, const QString& description)
{
    if (filtered_atom_tags.empty() && filtered_bond_tags.empty() &&
        filtered_non_molecular_tags.empty() && filtered_s_group_tags.empty()) {
        // nothing to select or deselect
        return;
    }

    auto redo = [this, filtered_atom_tags, filtered_bond_tags,
                 filtered_s_group_tags, filtered_non_molecular_tags,
                 to_select]() {
        setSelectionCommandFunc(filtered_atom_tags, filtered_bond_tags,
                                filtered_s_group_tags,
                                filtered_non_molecular_tags, to_select);
    };
    auto undo = [this, filtered_atom_tags, filtered_bond_tags,
                 filtered_s_group_tags, filtered_non_molecular_tags,
                 to_select]() {
        setSelectionCommandFunc(filtered_atom_tags, filtered_bond_tags,
                                filtered_s_group_tags,
                                filtered_non_molecular_tags, !to_select);
    };
    doCommand(redo, undo, description);
}

void MolModel::clearSelection()
{
    if (!hasSelection()) {
        // nothing to clear
        return;
    }
    auto sel_atom_tags = m_selected_atom_tags;
    auto sel_bond_tags = m_selected_bond_tags;
    auto sel_s_group_tags = m_selected_s_group_tags;
    auto sel_non_molecular_tags = m_selected_non_molecular_tags;

    QString desc = QString("Clear selection");
    auto redo = [this]() { clearSelectionCommandFunc(); };
    auto undo = [this, sel_atom_tags, sel_bond_tags, sel_s_group_tags,
                 sel_non_molecular_tags]() {
        m_selected_atom_tags = sel_atom_tags;
        m_selected_bond_tags = sel_bond_tags;
        m_selected_s_group_tags = sel_s_group_tags;
        m_selected_non_molecular_tags = sel_non_molecular_tags;
    };
    doCommand(redo, undo, desc);
}

void MolModel::selectAll()
{
    auto [atom_tags_to_select, bond_tags_to_select, s_group_tags_to_select,
          non_molecular_tags_to_select] = getAllUnselectedTags();
    doSelectionCommand(atom_tags_to_select, bond_tags_to_select,
                       s_group_tags_to_select, non_molecular_tags_to_select,
                       true, "Select all");
}

std::tuple<std::unordered_set<AtomTag>, std::unordered_set<BondTag>,
           std::unordered_set<SGroupTag>, std::unordered_set<NonMolecularTag>>
MolModel::getAllUnselectedTags()
{
    std::unordered_set<AtomTag> deselected_atoms;
    for (const auto& [atom_tag_int, atom] : *m_mol.getAtomBookmarks()) {
        auto atom_tag = AtomTag(atom_tag_int);
        if (m_selected_atom_tags.find(atom_tag) == m_selected_atom_tags.end()) {
            deselected_atoms.insert(atom_tag);
        }
    }
    std::unordered_set<BondTag> deselected_bonds;
    for (const auto& [bond_tag_int, bond] : *m_mol.getBondBookmarks()) {
        auto bond_tag = BondTag(bond_tag_int);
        if (m_selected_bond_tags.find(bond_tag) == m_selected_bond_tags.end()) {
            deselected_bonds.insert(bond_tag);
        }
    }
    std::unordered_set<SGroupTag> deselected_s_groups;
    for (const auto& s_group : getSubstanceGroups(m_mol)) {
        auto s_group_tag = getTagForSGroup(s_group);
        if (m_selected_s_group_tags.find(s_group_tag) ==
            m_selected_s_group_tags.end()) {
            deselected_s_groups.insert(s_group_tag);
        }
    }
    std::unordered_set<NonMolecularTag> deselected_non_molecular_objects;
    for (const auto* non_molecular_obj : getNonMolecularObjects()) {
        auto tag = non_molecular_obj->getTag();
        if (m_selected_non_molecular_tags.find(tag) ==
            m_selected_non_molecular_tags.end()) {
            deselected_non_molecular_objects.insert(tag);
        }
    }
    return {deselected_atoms, deselected_bonds, deselected_s_groups,
            deselected_non_molecular_objects};
}

void MolModel::invertSelection()
{
    auto [atom_tags_to_select, bond_tags_to_select, s_group_tags_to_select,
          non_molecular_tags_to_select] = getAllUnselectedTags();
    auto undo_macro_raii = createUndoMacro("Invert selection");
    doSelectionCommand(m_selected_atom_tags, m_selected_bond_tags,
                       m_selected_s_group_tags, m_selected_non_molecular_tags,
                       false, "Deselect");
    doSelectionCommand(atom_tags_to_select, bond_tags_to_select,
                       s_group_tags_to_select, non_molecular_tags_to_select,
                       true, "Select");
}

void MolModel::removeSelected()
{
    remove(getSelectedAtoms(), getSelectedBonds(), getSelectedSGroups(),
           getSelectedNonMolecularObjects());
}

void MolModel::mutateAtoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                           const Element element)
{
    unsigned int atomic_num = static_cast<unsigned int>(element);
    auto create_atom = std::bind(make_new_atom, atomic_num);
    mutateAtoms(atoms, create_atom);
}

void MolModel::mutateAtoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                           const AtomQuery atom_query)
{
    auto query_func = ATOM_TOOL_QUERY_MAP.at(atom_query);
    auto query_func_ptr =
        std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY>(query_func());
    auto create_atom = std::bind(make_new_query_atom, query_func_ptr);
    mutateAtoms(atoms, create_atom);
}

void MolModel::mutateAtoms(
    const std::unordered_set<const RDKit::Atom*>& from_atoms,
    const RDKit::Atom& to_atom)
{
    auto create_atom = [to_atom]() {
        return std::make_shared<RDKit::Atom>(to_atom);
    };
    mutateAtoms(from_atoms, create_atom);
}

void MolModel::mutateAtoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                           const AtomFunc& create_atom)
{
    if (atoms.empty()) {
        return;
    }
    auto cmd_func = [this, atoms, create_atom]() {
        for (auto atom : atoms) {
            mutateAtomCommandFunc(getTagForAtom(atom), create_atom);
        }
    };
    doCommandUsingSnapshots(cmd_func, "Mutate atoms", WhatChanged::MOLECULE);
}

void MolModel::mutateRGroups(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    auto next_r_group_num = get_next_r_group_numbers(&m_mol, 1)[0];
    mutateRGroups(atoms, next_r_group_num);
}

void MolModel::mutateRGroups(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const unsigned int r_group_num)
{
    auto cmd_func = [this, atoms, r_group_num]() {
        auto create_atom =
            std::bind(rdkit_extensions::make_new_r_group, r_group_num);
        for (auto atom : atoms) {
            mutateAtomCommandFunc(getTagForAtom(atom), create_atom);
        }
    };
    doCommandUsingSnapshots(cmd_func, "Mutate atoms", WhatChanged::MOLECULE);
}

void MolModel::mutateBonds(const std::unordered_set<const RDKit::Bond*>& bonds,
                           BondTool bond_tool)
{
    if (bonds.empty()) {
        return;
    }

    if (bond_tool == BondTool::ATOM_CHAIN) {
        return; // nothing to do for the atom chain tool
    }

    if (BOND_TOOL_BOND_MAP.count(bond_tool)) {
        auto [bond_type, bond_dir, flip_matching_bonds, cursor_hint_path] =
            BOND_TOOL_BOND_MAP.at(bond_tool);
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        mutateBonds(bonds, create_bond, flip_matching_bonds);

    } else if (BOND_TOOL_QUERY_MAP.count(bond_tool)) {
        auto query_func = BOND_TOOL_QUERY_MAP.at(bond_tool);
        auto bond_query =
            std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY>(query_func());
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        mutateBonds(bonds, create_bond);

    } else {
        throw std::runtime_error(
            fmt::format("Cannot mutate bonds for the given type {}",
                        static_cast<int>(bond_tool)));
    }
}

void MolModel::mutateBonds(const std::unordered_set<const RDKit::Bond*>& bonds,
                           const BondFunc& create_bond,
                           bool flip_matching_bonds)
{
    if (bonds.empty()) {
        return;
    }

    auto new_bond = create_bond();
    std::unordered_set<const RDKit::Bond*> matching_bonds;
    for (auto bond : bonds) {
        if (bond->getBondType() == new_bond->getBondType() &&
            bond->getBondDir() == new_bond->getBondDir()) {
            matching_bonds.insert(bond);
        }
    }

    auto cmd_func = [this, bonds, create_bond, flip_matching_bonds,
                     matching_bonds]() {
        for (auto bond : bonds) {
            mutateBondCommandFunc(getTagForBond(bond), create_bond);
        }
        if (flip_matching_bonds) {
            for (auto bond : matching_bonds) {
                flipBondCommandFunc(bond);
            }
        }
    };
    doCommandUsingSnapshots(cmd_func, "Mutate bonds", WhatChanged::MOLECULE);
}

void MolModel::adjustChargeOnAtoms(
    const std::unordered_set<const RDKit::Atom*>& atoms, const int increment_by)
{
    if (atoms.empty()) {
        return;
    }
    auto cmd_func = [this, increment_by, atoms]() {
        for (auto atom : atoms) {
            getMutableAtom(atom)->setFormalCharge(atom->getFormalCharge() +
                                                  increment_by);
        }
    };
    doCommandUsingSnapshots(cmd_func, "Set atomic charge",
                            WhatChanged::MOLECULE);
}

void MolModel::adjustRadicalElectronsOnAtoms(
    const std::unordered_set<const RDKit::Atom*>& atoms, int increment_by)
{
    if (atoms.empty()) {
        return;
    }
    auto cmd_func = [this, atoms, increment_by]() {
        for (auto atom : atoms) {
            getMutableAtom(atom)->setNumRadicalElectrons(
                atom->getNumRadicalElectrons() + increment_by);
        }
    };
    doCommandUsingSnapshots(cmd_func, "Set unpaired electrons",
                            WhatChanged::MOLECULE);
}

void MolModel::toggleExplicitHsOnAtoms(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    if (atoms.empty()) {
        return;
    }
    updateExplicitHs(ExplicitHActions::TOGGLE, atoms);
}

void MolModel::setTagForAtom(RDKit::Atom* const atom, const AtomTag atom_tag)
{
    m_mol.setAtomBookmark(atom, atom_tag);
    atom->setProp<int>(TAG_PROPERTY, atom_tag);
}

AtomTag MolModel::getTagForAtom(const RDKit::Atom* const atom,
                                const bool allow_null)
{
    if (atom == nullptr) {
        if (allow_null) {
            return AtomTag(-1);
        }
        throw std::runtime_error("Cannot pass nullptr to getTagForAtom "
                                 "unless allow_null is true");
    }
    return AtomTag(atom->getProp<int>(TAG_PROPERTY));
}

void MolModel::setTagForBond(RDKit::Bond* const bond, const BondTag bond_tag)
{
    m_mol.setBondBookmark(bond, bond_tag);
    bond->setProp<int>(TAG_PROPERTY, bond_tag);
}

BondTag MolModel::getTagForBond(const RDKit::Bond* const bond)
{
    return BondTag(bond->getProp<int>(TAG_PROPERTY));
}

const RDKit::Atom* MolModel::getAtomFromTag(AtomTag atom_tag) const
{
    // RDKit is missing const versions of bookmark getters, even though it
    // has const atom getters that take an atom index.  To get around this,
    // we use const_cast.  (See SHARED-9673.)
    return const_cast<RDKit::RWMol*>(&m_mol)->getUniqueAtomWithBookmark(
        atom_tag);
}

const RDKit::Bond* MolModel::getBondFromTag(BondTag bond_tag) const
{
    // RDKit is missing const versions of bookmark getters, even though it
    // has const bond getters that take a bond index.  To get around this,
    // we use const_cast.  (See SHARED-9673.)
    return const_cast<RDKit::RWMol*>(&m_mol)->getUniqueBondWithBookmark(
        bond_tag);
}

SGroupTag MolModel::getTagForSGroup(const RDKit::SubstanceGroup& s_group) const
{
    return SGroupTag(s_group.getProp<int>(TAG_PROPERTY));
}

void MolModel::setTagForSGroup(const RDKit::SubstanceGroup& s_group,
                               const SGroupTag s_group_tag)
{
    s_group.setProp<int>(TAG_PROPERTY, s_group_tag);
}

RDKit::SubstanceGroup
MolModel::getSGroupFromTag(const SGroupTag s_group_tag) const
{
    for (const auto& s_group : getSubstanceGroups(m_mol)) {
        if (getTagForSGroup(s_group) == s_group_tag) {
            return s_group;
        }
    }
    throw std::out_of_range("No S-group found");
}

RDKit::Atom* MolModel::getMutableAtom(const RDKit::Atom* const atom)
{
    return m_mol.getUniqueAtomWithBookmark(getTagForAtom(atom));
}
RDKit::Bond* MolModel::getMutableBond(const RDKit::Bond* const bond)
{
    return m_mol.getUniqueBondWithBookmark(getTagForBond(bond));
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

void MolModel::addAtomChainCommandFunc(
    const AtomFunc create_atom, const std::vector<RDGeom::Point3D>& coords,
    const BondFunc create_bond, const AtomTag bound_to_atom_tag)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Atom* prev_atom = nullptr;
    if (bound_to_atom_tag >= 0) {
        prev_atom = m_mol.getUniqueAtomWithBookmark(bound_to_atom_tag);
    }

    for (auto cur_coords : coords) {
        auto atom_shared = create_atom();
        auto* atom = atom_shared.get();
        unsigned int atom_index =
            m_mol.addAtom(atom, /* updateLabel = */ false);
        atom = m_mol.getAtomWithIdx(atom_index);
        setTagForAtom(atom, m_next_atom_tag++);
        m_mol.getConformer().setAtomPos(atom_index, cur_coords);

        if (prev_atom != nullptr) {
            auto bond_shared = create_bond();
            auto* bond = bond_shared.get();
            bond->setOwningMol(m_mol);
            bond->setBeginAtom(prev_atom);
            bond->setEndAtom(atom);
            unsigned int bond_index = m_mol.addBond(bond) - 1;
            bond = m_mol.getBondWithIdx(bond_index);
            setTagForBond(bond, m_next_bond_tag++);
        }
        prev_atom = atom;
    }
}

void MolModel::addNonMolecularObjectCommandFunc(const NonMolecularType& type,
                                                const RDGeom::Point3D& coords)
{
    Q_ASSERT(m_allow_edits);
    auto tag = m_next_non_molecular_tag++;
    if (type == NonMolecularType::RXN_ARROW) {
        m_arrow = NonMolecularObject(type, coords, tag);
    } else { // type == NonMolecularType::RXN_PLUS
        m_pluses.emplace_back(type, coords, tag);
    }
}

void MolModel::updateNonMolecularMetadata()
{
    m_tag_to_non_molecular_object.clear();
    for (auto& plus : m_pluses) {
        m_tag_to_non_molecular_object[plus.getTag()] = &plus;
    }
    if (m_arrow.has_value()) {
        m_tag_to_non_molecular_object[m_arrow->getTag()] = &m_arrow.value();
    }
}

bool MolModel::removeAtomCommandFunc(const AtomTag atom_tag)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
    bool selection_changed = m_selected_atom_tags.erase(atom_tag);
    // RDKit automatically deletes all bonds involving this atom, so we have
    // to remove those from the selection as well
    for (auto cur_bond : m_mol.atomBonds(atom)) {
        auto bond_tag = getTagForBond(cur_bond);
        if (m_selected_bond_tags.erase(bond_tag)) {
            selection_changed = true;
        }
    }
    m_mol.removeAtom(atom);
    return selection_changed;
}

void MolModel::addBondCommandFunc(const AtomTag start_atom_tag,
                                  const AtomTag end_atom_tag,
                                  const BondFunc create_bond)
{
    Q_ASSERT(m_allow_edits);
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
    setTagForBond(bond, m_next_bond_tag++);
}

bool MolModel::removeBondCommandFunc(const BondTag bond_tag,
                                     const AtomTag start_atom_tag,
                                     const AtomTag end_atom_tag)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Atom* start_atom = m_mol.getUniqueAtomWithBookmark(start_atom_tag);
    RDKit::Atom* end_atom = m_mol.getUniqueAtomWithBookmark(end_atom_tag);
    bool selection_changed = m_selected_bond_tags.erase(bond_tag);
    m_mol.removeBond(start_atom->getIdx(), end_atom->getIdx());
    return selection_changed;
}

bool MolModel::removeNonMolecularObjectCommandFunc(
    const NonMolecularTag cur_tag)
{
    auto* obj = m_tag_to_non_molecular_object.at(cur_tag);
    bool selected = m_selected_non_molecular_tags.erase(cur_tag);
    if (obj->getType() == NonMolecularType::RXN_ARROW) {
        m_arrow.reset();
    } else if (obj->getType() == NonMolecularType::RXN_PLUS) {
        auto iter_to_erase =
            std::find_if(m_pluses.begin(), m_pluses.end(),
                         [obj](const NonMolecularObject& cur_plus) {
                             return &cur_plus == obj;
                         });
        m_pluses.erase(iter_to_erase);
    }
    return selected;
}

void MolModel::removeCommandFunc(
    const std::vector<AtomTag>& atom_tags,
    const std::vector<std::tuple<BondTag, AtomTag, AtomTag>>&
        bond_tags_with_atoms,
    const std::unordered_set<const RDKit::SubstanceGroup*>& s_groups,
    const std::vector<NonMolecularTag>& non_molecular_tags)
{
    Q_ASSERT(m_allow_edits);
    // we have to determine whether we're deleting an attachment point before we
    // delete any of the bonds, since is_attachment_point() will return false
    // for unbound atoms
    bool attachment_point_deleted =
        std::any_of(atom_tags.begin(), atom_tags.end(), [this](AtomTag tag) {
            return is_attachment_point(getAtomFromTag(tag));
        });

    // handle the S-groups first so that they don't get implicitly deleted when
    // we remove an atom or a bond
    removeSGroupsCommandFunc(s_groups);
    deselectSGroupsThatWillBeImplicitlyDeleted(atom_tags, bond_tags_with_atoms);

    // then remove the bonds so that they don't get implicitly deleted when we
    // remove an atom
    for (auto [bond_tag, start_atom_tag, end_atom_tag] : bond_tags_with_atoms) {
        removeBondCommandFunc(bond_tag, start_atom_tag, end_atom_tag);
    }
    for (auto cur_atom_tag : atom_tags) {
        removeAtomCommandFunc(cur_atom_tag);
    }
    if (attachment_point_deleted) {
        renumber_attachment_points(&m_mol);
    }
    for (auto cur_tag : non_molecular_tags) {
        removeNonMolecularObjectCommandFunc(cur_tag);
    }
}

void MolModel::removeSGroupsCommandFunc(
    const std::unordered_set<const RDKit::SubstanceGroup*>& s_groups)
{
    Q_ASSERT(m_allow_edits);
    if (s_groups.empty()) {
        return;
    }
    for (auto* s_group : s_groups) {
        auto tag = getTagForSGroup(*s_group);
        m_selected_s_group_tags.erase(tag);
    }
    rdkit_extensions::remove_sgroups_from_molecule(m_mol, s_groups);
}

void MolModel::deselectSGroupsThatWillBeImplicitlyDeleted(
    const std::vector<AtomTag>& atom_tags,
    const std::vector<std::tuple<BondTag, AtomTag, AtomTag>>&
        bond_tags_with_atoms)
{
    Q_ASSERT(m_allow_edits);
    if (m_selected_s_group_tags.empty()) {
        // there are no selected S-groups, so we don't have to worry about
        // deselecting any of them
        return;
    }

    // build sets of indices for all of the atoms and bonds that are going to be
    // removed
    std::unordered_set<unsigned int> atom_idxs_to_delete;
    std::transform(
        atom_tags.begin(), atom_tags.end(),
        std::inserter(atom_idxs_to_delete, atom_idxs_to_delete.begin()),
        [this](AtomTag tag) { return getAtomFromTag(tag)->getIdx(); });
    std::unordered_set<unsigned int> bond_idxs_to_delete;
    std::transform(
        bond_tags_with_atoms.begin(), bond_tags_with_atoms.end(),
        std::inserter(bond_idxs_to_delete, bond_idxs_to_delete.begin()),
        [this](auto tags) {
            auto bond_tag = std::get<0>(tags);
            return getBondFromTag(bond_tag)->getIdx();
        });

    // figure out which of the selected S-groups are going to be implicitly
    // deleted
    std::unordered_set<SGroupTag> tags_to_deselect;
    for (auto sel_tag : m_selected_s_group_tags) {
        auto cur_s_group = getSGroupFromTag(sel_tag);
        bool s_group_contains_atom_to_be_deleted = std::any_of(
            cur_s_group.getAtoms().begin(), cur_s_group.getAtoms().end(),
            [&atom_idxs_to_delete](auto idx) {
                return atom_idxs_to_delete.count(idx);
            });
        bool s_group_contains_bond_to_be_deleted = std::any_of(
            cur_s_group.getBonds().begin(), cur_s_group.getBonds().end(),
            [&bond_idxs_to_delete](auto idx) {
                return bond_idxs_to_delete.count(idx);
            });

        if (s_group_contains_atom_to_be_deleted ||
            s_group_contains_bond_to_be_deleted) {
            // this S-group is going to be implicitly deleted
            tags_to_deselect.insert(sel_tag);
            continue;
        }
    }
    // finally, actually do the deselection
    for (auto tag : tags_to_deselect) {
        m_selected_s_group_tags.erase(tag);
    }
}

void MolModel::addMolCommandFunc(const RDKit::ROMol& mol)
{
    Q_ASSERT(m_allow_edits);
    // get the starting index for the atoms and bonds to be inserted
    unsigned int atom_index = m_mol.getNumAtoms();
    unsigned int bond_index = m_mol.getNumBonds();
    m_mol.insertMol(mol);
    unsigned int new_num_atoms = m_mol.getNumAtoms();
    unsigned int new_num_bonds = m_mol.getNumBonds();

    // insertMol only copies the atoms and bonds, so we need to add the sgroups
    for (auto sgroup : getSubstanceGroups(mol)) {
        // update the sgroup's atom and bond indices
        std::vector<unsigned int> new_atom_indices;
        std::transform(sgroup.getAtoms().begin(), sgroup.getAtoms().end(),
                       std::back_inserter(new_atom_indices),
                       [atom_index](unsigned int old_index) {
                           return old_index + atom_index;
                       });
        sgroup.setAtoms(new_atom_indices);
        std::vector<unsigned int> new_bond_indices;
        std::transform(sgroup.getBonds().begin(), sgroup.getBonds().end(),
                       std::back_inserter(new_bond_indices),
                       [bond_index](unsigned int old_index) {
                           return old_index + bond_index;
                       });
        sgroup.setBonds(new_bond_indices);
        setTagForSGroup(sgroup, m_next_s_group_tag++);
        addSubstanceGroup(m_mol, sgroup);
    }
    bool attachment_point_added = false;
    for (; atom_index < new_num_atoms; ++atom_index) {
        RDKit::Atom* atom = m_mol.getAtomWithIdx(atom_index);
        setTagForAtom(atom, m_next_atom_tag++);
        attachment_point_added =
            attachment_point_added || is_attachment_point(atom);
    }

    for (; bond_index < new_num_bonds; ++bond_index) {
        RDKit::Bond* bond = m_mol.getBondWithIdx(bond_index);
        setTagForBond(bond, m_next_bond_tag++);
    }
    if (attachment_point_added) {
        renumber_attachment_points(&m_mol);
    }
}

void MolModel::addReactionCommandFunc(const RDKit::ChemicalReaction& reaction)
{
    unsigned int num_reactants = reaction.getNumReactantTemplates();
    unsigned int num_products = reaction.getNumProductTemplates();
    std::vector<RDKit::ROMOL_SPTR> mols_to_add;
    mols_to_add.reserve(num_reactants + num_products);
    auto reactants_it = reaction.getReactants();
    mols_to_add.insert(mols_to_add.end(), reactants_it.begin(),
                       reactants_it.end());
    auto products_it = reaction.getProducts();
    mols_to_add.insert(mols_to_add.end(), products_it.begin(),
                       products_it.end());

    // insert the first molecule without drawing anything before it
    auto first_reaction_mol = mols_to_add[0];
    if (!m_mol.getNumAtoms()) {
        center_on_origin(*first_reaction_mol);
    } else {
        move_molecule_to_the_right_of(*first_reaction_mol, m_mol,
                                      IMPORT_SPACING);
    }
    addMolCommandFunc(*first_reaction_mol);
    // insert the rest of the molecules and draw either a plus sign or an arrow
    // before each one
    for (size_t i = 1; i < mols_to_add.size(); ++i) {
        auto reaction_mol = mols_to_add[i];
        auto prev_mol = mols_to_add[i - 1];
        double spacing;
        NonMolecularType non_mol_type;
        if (i == num_reactants) {
            // this is the first product, so we draw an arrow before it
            spacing = ARROW_SPACING;
            non_mol_type = NonMolecularType::RXN_ARROW;
        } else {
            // otherwise, we draw a plus sign before it
            spacing = PLUS_SPACING;
            non_mol_type = NonMolecularType::RXN_PLUS;
        }
        auto midpoint =
            move_molecule_to_the_right_of(*reaction_mol, *prev_mol, spacing);
        addNonMolecularObjectCommandFunc(non_mol_type, midpoint);
        addMolCommandFunc(*reaction_mol);
    }
}

void MolModel::mutateAtomCommandFunc(const AtomTag atom_tag,
                                     const AtomFunc create_atom)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
    int atom_index = atom->getIdx();
    auto new_atom = create_atom();
    m_mol.replaceAtom(atom_index, new_atom.get());
    // replaceAtom creates a copy, so we need to fetch the "real" new atom
    auto* mutated_atom = m_mol.getAtomWithIdx(atom_index);
    // The bookmark is automatically updated, but the property is not (unless we
    // passed preserveProps = true, but that overwrites the R-group property)
    mutated_atom->setProp<int>(TAG_PROPERTY, atom_tag);
}

void MolModel::mutateBondCommandFunc(const BondTag bond_tag,
                                     const BondFunc create_bond)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Bond* bond = m_mol.getUniqueBondWithBookmark(bond_tag);
    int bond_index = bond->getIdx();
    auto new_bond = create_bond();
    m_mol.replaceBond(bond_index, new_bond.get());
    // replaceBond creates a copy, so we need to fetch the "real" new bond
    auto* mutated_bond = m_mol.getBondWithIdx(bond_index);
    // The bookmark is automatically updated, but we have to manually copy the
    // bond tag property
    mutated_bond->setProp<int>(TAG_PROPERTY, bond_tag);
}

void MolModel::flipBondCommandFunc(const RDKit::Bond* bond)
{
    Q_ASSERT(m_allow_edits);
    auto mol_bond = getMutableBond(bond);
    unsigned int orig_begin = mol_bond->getBeginAtomIdx();
    unsigned int orig_end = mol_bond->getEndAtomIdx();
    mol_bond->setEndAtomIdx(orig_begin);
    mol_bond->setBeginAtomIdx(orig_end);
}

void MolModel::setCoordinates(
    const std::vector<AtomTag>& atom_tags,
    const std::vector<RDGeom::Point3D>& atom_coords,
    const std::vector<NonMolecularTag>& non_mol_tags,
    const std::vector<RDGeom::Point3D>& non_mol_coords)
{
    Q_ASSERT(m_allow_edits);
    if (atom_tags.size() != atom_coords.size() ||
        non_mol_tags.size() != non_mol_coords.size()) {
        throw std::invalid_argument("setCoordinates: tags and coords must have "
                                    "the same size");
    }

    for (auto const& [cur_tag, cur_coords] :
         boost::combine(atom_tags, atom_coords)) {
        RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(cur_tag);
        m_mol.getConformer().setAtomPos(atom->getIdx(), cur_coords);
    }
    for (auto const& [cur_non_mol_tag, cur_coords] :
         boost::combine(non_mol_tags, non_mol_coords)) {
        auto* non_mol_obj = m_tag_to_non_molecular_object.at(cur_non_mol_tag);
        non_mol_obj->setCoords(cur_coords);
    }

    // update the brackets for the sgroups with the new coordinates; this is
    // otherwise called when WhatChanged::MOLECULE but needs to be explicitly
    // done here to update the bracket positions
    rdkit_extensions::update_s_group_brackets(m_mol);

    emit coordinatesChanged();
}

void MolModel::clearCommandFunc()
{
    Q_ASSERT(m_allow_edits);
    m_mol.clear();
    // clear() removes the conformer, so we have to reinitialize m_mol
    initializeMol();

    m_pluses.clear();
    m_arrow.reset();

    m_selected_atom_tags.clear();
    m_selected_bond_tags.clear();
    m_selected_non_molecular_tags.clear();
}

void MolModel::setSelectionCommandFunc(
    const std::unordered_set<AtomTag>& atom_tags,
    const std::unordered_set<BondTag>& bond_tags,
    const std::unordered_set<SGroupTag>& s_group_tags,
    const std::unordered_set<NonMolecularTag>& non_molecular_tags,
    const bool selected)
{
    Q_ASSERT(m_allow_edits);
    if (selected) {
        for (auto cur_atom_tag : atom_tags) {
            m_selected_atom_tags.insert(cur_atom_tag);
        }
        for (auto cur_bond_tag : bond_tags) {
            m_selected_bond_tags.insert(cur_bond_tag);
        }
        for (auto cur_tag : s_group_tags) {
            m_selected_s_group_tags.insert(cur_tag);
        }
        for (auto cur_tag : non_molecular_tags) {
            m_selected_non_molecular_tags.insert(cur_tag);
        }
    } else {
        for (auto cur_atom_tag : atom_tags) {
            m_selected_atom_tags.erase(cur_atom_tag);
        }
        for (auto cur_bond_tag : bond_tags) {
            m_selected_bond_tags.erase(cur_bond_tag);
        }
        for (auto cur_tag : s_group_tags) {
            m_selected_s_group_tags.erase(cur_tag);
        }
        for (auto cur_tag : non_molecular_tags) {
            m_selected_non_molecular_tags.erase(cur_tag);
        }
    }
    emit selectionChanged();
}

void MolModel::clearSelectionCommandFunc()
{
    Q_ASSERT(m_allow_edits);
    bool selection_changed = hasSelection();
    m_selected_atom_tags.clear();
    m_selected_bond_tags.clear();
    m_selected_s_group_tags.clear();
    m_selected_non_molecular_tags.clear();
    if (selection_changed) {
        emit selectionChanged();
    }
}

static std::vector<unsigned int>
get_atom_indices(const std::unordered_set<const RDKit::Atom*>& atoms)
{
    std::vector<unsigned int> indices;
    std::transform(atoms.begin(), atoms.end(), std::back_inserter(indices),
                   [](auto atom) { return atom->getIdx(); });
    return indices;
}

void MolModel::addExplicitHsCommandFunc(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    Q_ASSERT(m_allow_edits);
    rdkit_extensions::addHs(m_mol, get_atom_indices(atoms));

    // Add tags for new atoms/bonds so selection continues to work in MolModel
    for (auto atom : m_mol.atoms()) {
        if (!atom->hasProp(TAG_PROPERTY)) {
            setTagForAtom(atom, m_next_atom_tag++);
            auto bond = *(m_mol.atomBonds(atom).begin());
            setTagForBond(bond, m_next_bond_tag++);
        }
    }
}

void MolModel::removeExplicitHsCommandFunc(
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    Q_ASSERT(m_allow_edits);
    if (atoms.empty()) {
        rdkit_extensions::removeHs(m_mol);
    } else {
        rdkit_extensions::removeHs(m_mol, get_atom_indices(atoms));
    }

    // deselect any atoms and bonds that were removed
    std::unordered_set<AtomTag> new_atom_selection;
    for (auto atom_tag : m_selected_atom_tags) {
        if (m_mol.hasAtomBookmark(atom_tag)) {
            new_atom_selection.insert(atom_tag);
        }
    }
    m_selected_atom_tags = new_atom_selection;

    std::unordered_set<BondTag> new_bond_selection;
    for (auto bond_tag : m_selected_bond_tags) {
        if (m_mol.hasBondBookmark(bond_tag)) {
            new_bond_selection.insert(bond_tag);
        }
    }
    m_selected_bond_tags = new_bond_selection;
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/model/mol_model.moc"
