#include "schrodinger/sketcher/model/mol_model.h"

#include <algorithm>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
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
#include "schrodinger/rdkit_extensions/coord_utils.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"
#include "schrodinger/sketcher/molviewer/constants.h"
#include "schrodinger/sketcher/molviewer/coord_utils.h"
#include "schrodinger/sketcher/rdkit/fragment.h"
#include "schrodinger/sketcher/rdkit/molops.h"
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

MolModelSnapshot::MolModelSnapshot(
    const RDKit::RWMol& mol, const std::vector<NonMolecularObject>& pluses,
    const std::optional<NonMolecularObject>& arrow,
    const std::unordered_set<int>& selected_atom_tags,
    const std::unordered_set<int>& selected_bond_tags,
    const std::unordered_set<int>& selected_non_molecular_tags) :
    m_mol(mol),
    m_pluses(pluses),
    m_arrow(arrow),
    m_selected_atom_tags(selected_atom_tags),
    m_selected_bond_tags(selected_bond_tags),
    m_selected_non_molecular_tags(selected_non_molecular_tags)
{
}

bool MolModelSnapshot::isSelectionIdentical(const MolModelSnapshot& other)
{
    return m_selected_atom_tags == other.m_selected_atom_tags &&
           m_selected_bond_tags == other.m_selected_bond_tags &&
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
    // RDKit takes ownership of this conformer, so we don't have to worry about
    // deleting it
    m_mol.addConformer(conf, true);
}

const RDKit::ROMol* MolModel::getMol() const
{
    // TODO: add API to get selected mol
    // TODO: add API to export selection as atom/bond properties
    return &m_mol;
}

std::string MolModel::getMolText(const Format format) const
{
    return rdkit_extensions::to_string(*getMol(), format);
}

std::shared_ptr<RDKit::ChemicalReaction> MolModel::getReaction() const
{
    if (!m_arrow.has_value()) {
        throw std::runtime_error("No reaction arrow found.");
    }
    // get the midpoint of the arrow
    auto arrow_x = m_arrow->getCoords().x;
    auto all_mols =
        RDKit::MolOps::getMolFrags(m_mol, /* sanitizeFrags = */ false);
    auto reaction = std::make_shared<RDKit::ChemicalReaction>();
    for (auto cur_mol : all_mols) {
        auto cur_centroid = find_centroid(*cur_mol);
        if (cur_centroid.x <= arrow_x) {
            reaction->addReactantTemplate(cur_mol);
        } else {
            reaction->addProductTemplate(cur_mol);
        }
    }
    if (reaction->getReactants().empty() || reaction->getProducts().empty()) {
        throw std::runtime_error("Incomplete reactions cannot be copied.");
    }
    return reaction;
}

std::string MolModel::getReactionText(const Format format) const
{
    return rdkit_extensions::to_string(*getReaction(), format);
}

bool MolModel::isEmpty() const
{
    return !m_mol.getNumAtoms() && m_pluses.empty() && !m_arrow.has_value();
}

bool MolModel::hasSelection() const
{
    return !(m_selected_atom_tags.empty() && m_selected_bond_tags.empty() &&
             m_selected_non_molecular_tags.empty());
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

std::unordered_set<const NonMolecularObject*>
MolModel::getSelectedNonMolecularObjects() const
{
    std::unordered_set<const NonMolecularObject*> selected;

    for (int tag : m_selected_non_molecular_tags) {
        selected.insert(m_tag_to_non_molecular_object.at(tag));
    }
    return selected;
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
        update_molecule_metadata(m_mol);
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
                            m_selected_bond_tags,
                            m_selected_non_molecular_tags);
}

void MolModel::restoreSnapshot(const MolModelSnapshot& snapshot,
                               const WhatChangedType what_changed,
                               const bool selection_changed,
                               const bool arrow_added)
{
    Q_ASSERT(m_allow_edits);
    if (what_changed & WhatChanged::MOLECULE) {
        m_mol = snapshot.m_mol;
        // we don't need to call update_molecule_metadata since all metadata was
        // updated before we took the snapshot
    }
    if (what_changed & WhatChanged::NON_MOL_OBJS) {
        m_pluses = snapshot.m_pluses;
        m_arrow = snapshot.m_arrow;
        updateNonMolecularMetadata();
    }
    if (selection_changed) {
        m_selected_atom_tags = snapshot.m_selected_atom_tags;
        m_selected_bond_tags = snapshot.m_selected_bond_tags;
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
    int bound_to_atom_tag = getTagForAtom(bound_to_atom);
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
    int bound_to_atom_tag =
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
    int bound_to_atom_tag =
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
    int bound_to_atom_tag =
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
    int bound_to_atom_tag =
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
    int bound_to_atom_tag =
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
    int start_atom_tag = getTagForAtom(start_atom);
    int end_atom_tag = getTagForAtom(end_atom);

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
    int start_atom_tag = getTagForAtom(start_atom);
    int end_atom_tag = getTagForAtom(end_atom);

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
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects)
{
    WhatChangedType to_be_changed = WhatChanged::NOTHING;
    if (!(atoms.empty() && bonds.empty())) {
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

    std::vector<int> non_molecular_tags;
    for (auto* cur_obj : non_molecular_objects) {
        non_molecular_tags.push_back(cur_obj->getTag());
    }

    QString desc = QString("Erase");
    auto cmd_func = [this, atom_tags, bond_tags, non_molecular_tags]() {
        removeCommandFunc(atom_tags, bond_tags, non_molecular_tags);
    };
    doCommandUsingSnapshots(cmd_func, desc, to_be_changed);
}

void MolModel::addMol(RDKit::ROMol mol, const QString& description,
                      const bool reposition_mol)
{
    if (mol.getNumAtoms() == 0) {
        return;
    }
    if (reposition_mol) {
        // make sure that the molecule has coordinates
        rdkit_extensions::update_2d_coordinates(mol);
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

void MolModel::importFromText(const std::string& text, Format format)
{
    auto mol_or_reaction = text_to_mol_or_reaction(text, format);
    if (auto* reaction =
            std::get_if<boost::shared_ptr<RDKit::ChemicalReaction>>(
                &mol_or_reaction)) {
        addReaction(**reaction);
    } else {
        auto mol = std::get<boost::shared_ptr<RDKit::RWMol>>(mol_or_reaction);
        addMol(*mol);
    }
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
    std::vector<std::pair<unsigned int, AtomFunc>> mutations_to_core_atoms,
    std::vector<std::pair<unsigned int, BondFunc>> mutations_to_core_bonds,
    std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>
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

void MolModel::mutateAtom(const RDKit::Atom* const atom, const Element& element)
{
    int atom_tag = getTagForAtom(atom);
    unsigned int atomic_num = static_cast<unsigned int>(element);
    auto cmd_func = [this, atom_tag, atomic_num]() {
        auto create_atom = std::bind(make_new_atom, atomic_num);
        mutateAtomCommandFunc(atom_tag, create_atom);
    };
    doCommandUsingSnapshots(cmd_func, "Mutate atom", WhatChanged::MOLECULE);
}

void MolModel::setAtomCharge(const RDKit::Atom* const atom, int charge)
{
    int atom_tag = getTagForAtom(atom);
    auto cmd_func = [this, atom_tag, charge]() {
        setAtomChargeCommandFunc(atom_tag, charge);
    };
    doCommandUsingSnapshots(cmd_func, "Set atom charge", WhatChanged::MOLECULE);
}

void MolModel::setAtomMapping(
    const std::unordered_set<const RDKit::Atom*>& atoms, int mapping)
{
    if (atoms.empty()) {
        return;
    }
    std::unordered_set<int> atom_tags(atoms.size());
    for (const RDKit::Atom* atom : atoms) {
        atom_tags.insert(getTagForAtom(atom));
    }

    auto redo = [this, atom_tags, mapping]() {
        setAtomMappingCommandFunc(atom_tags, mapping);
    };

    doCommandUsingSnapshots(redo, "Set mapping number", WhatChanged::MOLECULE);
}

void MolModel::mutateAtom(
    const RDKit::Atom* const atom,
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query)
{
    int atom_tag = getTagForAtom(atom);
    auto cmd_func = [this, atom_tag, atom_query]() {
        auto create_atom = std::bind(make_new_query_atom, atom_query);
        mutateAtomCommandFunc(atom_tag, create_atom);
    };
    doCommandUsingSnapshots(cmd_func, "Mutate atom", WhatChanged::MOLECULE);
}

bool MolModel::hasAnyImplicitHs(
    const std::unordered_set<const RDKit::Atom*>& atoms) const
{
    // RDKit calls some hydrogens "explicit" but they are not really
    // explicit in the sense of being in the graph and we need to count them
    // in. We count all hydrogens that are not in the graph as implicit.
    return std::any_of(atoms.begin(), atoms.end(),
                       [](auto atom) { return atom->getTotalNumHs() > 0; });
}

void MolModel::updateExplicitHs(ExplicitHActions action,
                                std::unordered_set<const RDKit::Atom*> atoms)
{
    bool add_hs;
    if (action == ExplicitHActions::TOGGLE) {
        std::unordered_set<const RDKit::Atom*> atoms_to_check;
        if (atoms.empty()) {
            atoms_to_check = atoms;
        } else {
            auto all_atoms = m_mol.atoms();
            atoms_to_check.insert(all_atoms.begin(), all_atoms.end());
        }
        add_hs = hasAnyImplicitHs(atoms_to_check);
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

void MolModel::mutateRGroup(const RDKit::Atom* const atom,
                            const unsigned int r_group_num)
{
    int atom_tag = getTagForAtom(atom);
    auto cmd_func = [this, atom_tag, r_group_num]() {
        auto create_atom =
            std::bind(rdkit_extensions::make_new_r_group, r_group_num);
        mutateAtomCommandFunc(atom_tag, create_atom);
    };
    doCommandUsingSnapshots(cmd_func, "Mutate atom", WhatChanged::MOLECULE);
}

void MolModel::mutateBond(const RDKit::Bond* const bond,
                          const RDKit::Bond::BondType& bond_type,
                          const RDKit::Bond::BondDir& bond_dir)
{
    int bond_tag = getTagForBond(bond);
    auto cmd_func = [this, bond_tag, bond_type, bond_dir]() {
        auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
        mutateBondCommandFunc(bond_tag, create_bond);
    };
    doCommandUsingSnapshots(cmd_func, "Mutate bond", WhatChanged::MOLECULE);
}

void MolModel::mutateBond(
    const RDKit::Bond* const bond,
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    int bond_tag = getTagForBond(bond);
    auto cmd_func = [this, bond_tag, bond_query]() {
        auto create_bond = std::bind(make_new_query_bond, bond_query);
        mutateBondCommandFunc(bond_tag, create_bond);
    };
    doCommandUsingSnapshots(cmd_func, "Mutate bond", WhatChanged::MOLECULE);
}

void MolModel::flipBond(const RDKit::Bond* const bond)
{
    int bond_tag = getTagForBond(bond);
    auto cmd_func = [this, bond_tag]() { flipBondCommandFunc(bond_tag); };
    doCommandUsingSnapshots(cmd_func, "Flip bond", WhatChanged::MOLECULE);
}

void MolModel::regenerateCoordinates()
{
    auto cmd = [this]() {
        rdkit_extensions::compute2DCoords(m_mol);
        assign_CIP_labels(m_mol);
        emit moleculeChanged();
    };
    doCommandUsingSnapshots(cmd, "Clean Up Coordinates", WhatChanged::MOLECULE);
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

    // get tags
    std::vector<int> atom_tags(atom_vec.size());
    transform(atom_vec.begin(), atom_vec.end(), atom_tags.begin(),
              [this](const RDKit::Atom* atom) { return getTagForAtom(atom); });
    std::vector<int> non_mol_tags(non_mol_vec.size());
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
        typedef std::tuple<std::vector<int>, std::vector<RDGeom::Point3D>,
                           std::vector<RDGeom::Point3D>, std::vector<int>,
                           std::vector<RDGeom::Point3D>,
                           std::vector<RDGeom::Point3D>>
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

void MolModel::flipAroundSegment(
    const RDGeom::Point3D& p1, const RDGeom::Point3D& p2,
    const std::unordered_set<const RDKit::Atom*>& atoms)
{
    auto flip = [p1, p2](auto& coord) { coord = flip_point(coord, p1, p2); };
    transformCoordinatesWithFunction("Flip selection", flip, MergeId::NO_MERGE,
                                     atoms);
}

void MolModel::flipAllHorizontal()
{
    auto center = find_centroid(m_mol, getNonMolecularObjects());
    auto flip_x = [center](auto& coord) { coord.x = 2 * center.x - coord.x; };
    transformCoordinatesWithFunction("Flip all horizontal", flip_x);
}

void MolModel::flipAllVertical()
{
    auto center = find_centroid(m_mol, getNonMolecularObjects());
    auto flip_y = [center](auto& coord) { coord.y = 2 * center.y - coord.y; };
    transformCoordinatesWithFunction("Flip all vertical", flip_y);
}

void MolModel::select(
    const std::unordered_set<const RDKit::Atom*>& atoms,
    const std::unordered_set<const RDKit::Bond*>& bonds,
    const std::unordered_set<const NonMolecularObject*>& non_molecular_objects,
    const SelectMode select_mode)
{
    auto [expanded_atoms, expanded_bonds] =
        ensureCompleteAttachmentPoints(atoms, bonds);
    std::unordered_set<int> atom_tags;
    for (auto* cur_atom : expanded_atoms) {
        atom_tags.insert(getTagForAtom(cur_atom));
    }
    std::unordered_set<int> bond_tags;
    for (auto* cur_bond : expanded_bonds) {
        bond_tags.insert(getTagForBond(cur_bond));
    }
    std::unordered_set<int> non_molecular_tags;
    for (auto* cur_non_molecular_obj : non_molecular_objects) {
        non_molecular_tags.insert(cur_non_molecular_obj->getTag());
    }
    selectTags(atom_tags, bond_tags, non_molecular_tags, select_mode);
}

void MolModel::selectTags(const std::unordered_set<int>& atom_tags,
                          const std::unordered_set<int>& bond_tags,
                          const std::unordered_set<int>& non_molecular_tags,
                          const SelectMode select_mode)
{
    // note that we do not ensure that attachment points are complete here.
    // That happens in select, not selectTags
    bool no_tags_specified =
        atom_tags.empty() && bond_tags.empty() && non_molecular_tags.empty();
    if (select_mode == SelectMode::SELECT_ONLY) {
        if (no_tags_specified) {
            clearSelection();
        } else {
            auto undo_macro_raii = createUndoMacro("Select only");
            clearSelection();
            doSelectionCommand(atom_tags, bond_tags, non_molecular_tags, true,
                               "Select");
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
    auto [selected_non_molecular_tags, deselected_non_molecular_tags] =
        divideBySelected(non_molecular_tags, m_selected_non_molecular_tags);
    if (select_mode == SelectMode::SELECT) {
        // if we passed atom_tags and bond_tags to doSelectionCommand
        // instead of deselected_atom_tags and deselected_bond_tags, then
        // undoing the command could deselect atoms and bonds that should've
        // remained selected.
        doSelectionCommand(deselected_atom_tags, deselected_bond_tags,
                           deselected_non_molecular_tags, true, "Select");
    } else if (select_mode == SelectMode::DESELECT) {
        // if we passed atom_tags and bond_tags to doSelectionCommand
        // instead of selected_atom_tags and selected_bond_tags, then
        // undoing the command could select atoms and bonds that should've
        // remained deselected.
        doSelectionCommand(selected_atom_tags, selected_bond_tags,
                           selected_non_molecular_tags, false, "Deselect");
    } else { // select_mode == SelectMode::TOGGLE
        auto undo_macro_raii = createUndoMacro("Toggle selection");
        doSelectionCommand(deselected_atom_tags, deselected_bond_tags,
                           deselected_non_molecular_tags, true, "Select");
        doSelectionCommand(selected_atom_tags, selected_bond_tags,
                           selected_non_molecular_tags, false, "Deselect");
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
    const std::unordered_set<int>& filtered_bond_tags,
    const std::unordered_set<int>& filtered_non_molecular_tags,
    const bool to_select, const QString& description)
{
    if (filtered_atom_tags.empty() && filtered_bond_tags.empty() &&
        filtered_non_molecular_tags.empty()) {
        // nothing to select or deselect
        return;
    }

    auto redo = [this, filtered_atom_tags, filtered_bond_tags,
                 filtered_non_molecular_tags, to_select]() {
        setSelectionCommandFunc(filtered_atom_tags, filtered_bond_tags,
                                filtered_non_molecular_tags, to_select);
    };
    auto undo = [this, filtered_atom_tags, filtered_bond_tags,
                 filtered_non_molecular_tags, to_select]() {
        setSelectionCommandFunc(filtered_atom_tags, filtered_bond_tags,
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
    std::unordered_set<int> sel_atom_tags = m_selected_atom_tags;
    std::unordered_set<int> sel_bond_tags = m_selected_bond_tags;
    std::unordered_set<int> sel_non_molecular_tags =
        m_selected_non_molecular_tags;

    QString desc = QString("Clear selection");
    auto redo = [this]() { clearSelectionCommandFunc(); };
    auto undo = [this, sel_atom_tags, sel_bond_tags, sel_non_molecular_tags]() {
        m_selected_atom_tags = sel_atom_tags;
        m_selected_bond_tags = sel_bond_tags;
        m_selected_non_molecular_tags = sel_non_molecular_tags;
    };
    doCommand(redo, undo, desc);
}

void MolModel::selectAll()
{
    auto [atom_tags_to_select, bond_tags_to_select,
          non_molecular_tags_to_select] = getAllUnselectedTags();
    doSelectionCommand(atom_tags_to_select, bond_tags_to_select,
                       non_molecular_tags_to_select, true, "Select all");
}

std::tuple<std::unordered_set<int>, std::unordered_set<int>,
           std::unordered_set<int>>
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
    std::unordered_set<int> deselected_non_molecular_objects;
    for (const auto* non_molecular_obj : getNonMolecularObjects()) {
        int tag = non_molecular_obj->getTag();
        if (m_selected_non_molecular_tags.find(tag) ==
            m_selected_non_molecular_tags.end()) {
            deselected_non_molecular_objects.insert(tag);
        }
    }
    return {deselected_atoms, deselected_bonds,
            deselected_non_molecular_objects};
}

void MolModel::invertSelection()
{
    auto [atom_tags_to_select, bond_tags_to_select,
          non_molecular_tags_to_select] = getAllUnselectedTags();
    auto undo_macro_raii = createUndoMacro("Invert selection");
    doSelectionCommand(m_selected_atom_tags, m_selected_bond_tags,
                       m_selected_non_molecular_tags, false, "Deselect");
    doSelectionCommand(atom_tags_to_select, bond_tags_to_select,
                       non_molecular_tags_to_select, true, "Select");
}

void MolModel::removeSelected()
{
    remove(getSelectedAtoms(), getSelectedBonds(),
           getSelectedNonMolecularObjects());
}

void MolModel::mutateSelectedAtoms(const Element& element)
{
    unsigned int atomic_num = static_cast<unsigned int>(element);
    auto create_atom = std::bind(make_new_atom, atomic_num);
    mutateSelectedAtoms(create_atom);
}

void MolModel::mutateSelectedAtoms(
    const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query)
{
    auto create_atom = std::bind(make_new_query_atom, atom_query);
    mutateSelectedAtoms(create_atom);
}

void MolModel::mutateSelectedAtoms(const AtomFunc& create_atom)
{
    if (m_selected_atom_tags.empty()) {
        return;
    }
    auto cmd_func = [this, create_atom]() {
        for (auto atom_tag : m_selected_atom_tags) {
            mutateAtomCommandFunc(atom_tag, create_atom);
        }
    };
    doCommandUsingSnapshots(cmd_func, "Mutate atoms", WhatChanged::MOLECULE);
}

void MolModel::mutateSelectedBonds(const RDKit::Bond::BondType& bond_type,
                                   const RDKit::Bond::BondDir& bond_dir,
                                   bool flip_matching_bonds)
{
    std::unordered_set<int> matching_bond_tags, non_matching_bond_tags;
    for (auto* bond : getSelectedBonds()) {
        auto bond_tag = getTagForBond(bond);
        if (bond->getBondType() == bond_type &&
            bond->getBondDir() == bond_dir) {
            matching_bond_tags.insert(bond_tag);
        } else {
            non_matching_bond_tags.insert(bond_tag);
        }
    }
    if (non_matching_bond_tags.empty() &&
        (!flip_matching_bonds || matching_bond_tags.empty())) {
        // nothing to do
        return;
    }

    auto create_bond = std::bind(make_new_bond, bond_type, bond_dir);
    auto cmd_func = [this, create_bond, flip_matching_bonds, matching_bond_tags,
                     non_matching_bond_tags]() {
        for (auto bond_tag : non_matching_bond_tags) {
            mutateBondCommandFunc(bond_tag, create_bond);
        }
        if (flip_matching_bonds) {
            for (auto bond_tag : matching_bond_tags) {
                flipBondCommandFunc(bond_tag);
            }
        }
    };
    doCommandUsingSnapshots(cmd_func, "Mutate bonds", WhatChanged::MOLECULE);
}

void MolModel::mutateSelectedBonds(
    const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query)
{
    if (m_selected_bond_tags.empty()) {
        return;
    }
    auto create_bond = std::bind(make_new_query_bond, bond_query);
    auto cmd_func = [this, create_bond]() {
        for (auto bond_tag : m_selected_bond_tags) {
            mutateBondCommandFunc(bond_tag, create_bond);
        }
    };
    doCommandUsingSnapshots(cmd_func, "Mutate bonds", WhatChanged::MOLECULE);
}

void MolModel::adjustChargeOnSelectedAtoms(int increment_by)
{
    if (m_selected_atom_tags.empty()) {
        return;
    }
    auto cmd_func = [this, increment_by]() {
        Q_ASSERT(m_allow_edits);
        for (auto atom_tag : m_selected_atom_tags) {
            auto* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
            atom->setFormalCharge(atom->getFormalCharge() + increment_by);
        }
    };
    doCommandUsingSnapshots(cmd_func, "Set atom charge", WhatChanged::MOLECULE);
}

void MolModel::toggleExplicitHsOnSelectedAtoms()
{
    updateExplicitHs(ExplicitHActions::TOGGLE, getSelectedAtoms());
}

void MolModel::setTagForAtom(RDKit::Atom* const atom, const int atom_tag)
{
    m_mol.setAtomBookmark(atom, atom_tag);
    atom->setProp(TAG_PROPERTY, atom_tag);
}

int MolModel::getTagForAtom(const RDKit::Atom* const atom,
                            const bool allow_null)
{
    if (atom == nullptr) {
        if (allow_null) {
            return -1;
        }
        throw std::runtime_error("Cannot pass nullptr to getTagForAtom "
                                 "unless allow_null is true");
    }
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

void MolModel::addAtomChainCommandFunc(
    const AtomFunc create_atom, const std::vector<RDGeom::Point3D>& coords,
    const BondFunc create_bond, const int bound_to_atom_tag)
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

bool MolModel::removeAtomCommandFunc(const int atom_tag)
{
    Q_ASSERT(m_allow_edits);
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

void MolModel::addBondCommandFunc(const int start_atom_tag,
                                  const int end_atom_tag,
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

bool MolModel::removeBondCommandFunc(const int bond_tag,
                                     const int start_atom_tag,
                                     const int end_atom_tag)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Atom* start_atom = m_mol.getUniqueAtomWithBookmark(start_atom_tag);
    RDKit::Atom* end_atom = m_mol.getUniqueAtomWithBookmark(end_atom_tag);
    bool selection_changed = m_selected_bond_tags.erase(bond_tag);
    m_mol.removeBond(start_atom->getIdx(), end_atom->getIdx());
    return selection_changed;
}

bool MolModel::removeNonMolecularObjectCommandFunc(const int cur_tag)
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
    const std::vector<int>& atom_tags,
    const std::vector<std::tuple<int, int, int>>& bond_tags_with_atoms,
    const std::vector<int>& non_molecular_tags)
{
    Q_ASSERT(m_allow_edits);
    // we have to determine whether we're deleting an attachment point
    // before we delete any of the bonds, since is_attachment_point() will
    // return false for unbound atoms
    bool attachment_point_deleted =
        std::any_of(atom_tags.begin(), atom_tags.end(), [this](int tag) {
            return is_attachment_point(getAtomFromTag(tag));
        });
    // remove the bonds first so that they don't get implicitly deleted when
    // we remove an atom
    for (auto [bond_tag, start_atom_tag, end_atom_tag] : bond_tags_with_atoms) {
        removeBondCommandFunc(bond_tag, start_atom_tag, end_atom_tag);
    }
    for (int cur_atom_tag : atom_tags) {
        removeAtomCommandFunc(cur_atom_tag);
    }
    if (attachment_point_deleted) {
        renumber_attachment_points(&m_mol);
    }
    for (int cur_tag : non_molecular_tags) {
        removeNonMolecularObjectCommandFunc(cur_tag);
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

void MolModel::mutateAtomCommandFunc(const int atom_tag,
                                     const AtomFunc create_atom)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
    int atom_index = atom->getIdx();
    auto new_atom = create_atom();
    m_mol.replaceAtom(atom_index, new_atom.get());
    // replaceAtom creates a copy, so we need to fetch the "real" new atom
    auto* mutated_atom = m_mol.getAtomWithIdx(atom_index);
    // The bookmark is automatically updated, but the property is not
    // (unless we passed preserveProbs = true, but that overwrites the
    // R-group property)
    mutated_atom->setProp(TAG_PROPERTY, atom_tag);
}

void MolModel::mutateBondCommandFunc(const int bond_tag,
                                     const BondFunc create_bond)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Bond* bond = m_mol.getUniqueBondWithBookmark(bond_tag);
    int bond_index = bond->getIdx();
    auto new_bond = create_bond();
    m_mol.replaceBond(bond_index, new_bond.get());
    // replaceBond creates a copy, so we need to fetch the "real" new bond
    auto* mutated_bond = m_mol.getBondWithIdx(bond_index);
    // The bookmark is automatically updated, but we have to manually copy
    // the bond tag property
    mutated_bond->setProp(TAG_PROPERTY, bond_tag);
}

void MolModel::setAtomMappingCommandFunc(
    const std::unordered_set<int>& atom_tags, const int atom_mapping)
{
    Q_ASSERT(m_allow_edits);
    for (auto atom_tag : atom_tags) {
        RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
        atom->setAtomMapNum(atom_mapping);
    }
}

void MolModel::setAtomChargeCommandFunc(const int atom_tag, const int charge)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
    atom->setFormalCharge(charge);
}

void MolModel::flipBondCommandFunc(const int bond_tag)
{
    Q_ASSERT(m_allow_edits);
    RDKit::Bond* bond = m_mol.getUniqueBondWithBookmark(bond_tag);
    unsigned int orig_begin = bond->getBeginAtomIdx();
    unsigned int orig_end = bond->getEndAtomIdx();
    bond->setEndAtomIdx(orig_begin);
    bond->setBeginAtomIdx(orig_end);
}

void MolModel::setCoordinates(
    const std::vector<int>& atom_tags,
    const std::vector<RDGeom::Point3D>& atom_coords,
    const std::vector<int>& non_mol_tags,
    const std::vector<RDGeom::Point3D>& non_mol_coords)
{
    Q_ASSERT(m_allow_edits);
    if (atom_tags.size() != atom_coords.size() ||
        non_mol_tags.size() != non_mol_coords.size()) {
        throw std::invalid_argument("setCoordinates: tags and coords must have "
                                    "the same size");
    }

    int cur_tag;
    RDGeom::Point3D cur_coords;
    BOOST_FOREACH (boost::tie(cur_tag, cur_coords),
                   boost::combine(atom_tags, atom_coords)) {
        RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(cur_tag);
        m_mol.getConformer().setAtomPos(atom->getIdx(), cur_coords);
    }
    BOOST_FOREACH (boost::tie(cur_tag, cur_coords),
                   boost::combine(non_mol_tags, non_mol_coords)) {
        auto* non_mol_obj = m_tag_to_non_molecular_object.at(cur_tag);
        non_mol_obj->setCoords(cur_coords);
    }

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
    const std::unordered_set<int>& atom_tags,
    const std::unordered_set<int>& bond_tags,
    const std::unordered_set<int>& non_molecular_tags, const bool selected)
{
    Q_ASSERT(m_allow_edits);
    if (selected) {
        for (int cur_atom_tag : atom_tags) {
            m_selected_atom_tags.insert(cur_atom_tag);
        }
        for (int cur_bond_tag : bond_tags) {
            m_selected_bond_tags.insert(cur_bond_tag);
        }
        for (int cur_tag : non_molecular_tags) {
            m_selected_non_molecular_tags.insert(cur_tag);
        }
    } else {
        for (int cur_atom_tag : atom_tags) {
            m_selected_atom_tags.erase(cur_atom_tag);
        }
        for (int cur_bond_tag : bond_tags) {
            m_selected_bond_tags.erase(cur_bond_tag);
        }
        for (int cur_tag : non_molecular_tags) {
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
    std::unordered_set<int> new_atom_selection;
    for (auto atom_tag : m_selected_atom_tags) {
        if (m_mol.hasAtomBookmark(atom_tag)) {
            new_atom_selection.insert(atom_tag);
        }
    }
    m_selected_atom_tags = new_atom_selection;

    std::unordered_set<int> new_bond_selection;
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
