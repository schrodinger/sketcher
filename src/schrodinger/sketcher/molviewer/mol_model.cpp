#include "schrodinger/sketcher/molviewer/mol_model.h"

#include <boost/range/combine.hpp>

#include <GraphMol/Atom.h>
#include <GraphMol/Bond.h>
#include <GraphMol/Conformer.h>
#include <GraphMol/RWMol.h>

#include <QtGlobal>
#include <QUndoStack>

namespace schrodinger
{
namespace sketcher
{

MolModel::MolModel(QUndoStack* const undo_stack) :
    AbstractUndoableModel(undo_stack)
{
    auto* conf = new RDKit::Conformer();
    // RDKit takes ownership of this conformer, so we don't have to worry about
    // deleting it
    m_mol.addConformer(conf, true);
}

const RDKit::ROMol* MolModel::getMol() const
{
    return &m_mol;
}

void MolModel::doCommandWithMolUndo(const std::function<void()> redo,
                                    const QString& description)
{
    RDKit::ROMol undo_mol = m_mol;
    auto undo = [this, undo_mol]() {
        m_mol = undo_mol;
        emit moleculeChanged();
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

void MolModel::removeAtom(const RDKit::Atom* const atom)
{
    int atom_tag = getTagForAtom(atom);
    std::string element = atom->getSymbol();
    QString desc = QString("Remove %1").arg(QString::fromStdString(element));
    auto redo = [this, atom_tag]() { removeAtomFromCommand(atom_tag); };
    doCommandWithMolUndo(redo, desc);
}

void MolModel::removeBond(const RDKit::Bond* const bond)
{
    int bond_tag = getTagForBond(bond);
    int start_atom_tag = getTagForAtom(bond->getBeginAtom());
    int end_atom_tag = getTagForAtom(bond->getEndAtom());

    QString desc = QString("Remove bond");
    auto redo = [this, bond_tag, start_atom_tag, end_atom_tag]() {
        removeBondFromCommand(bond_tag, start_atom_tag, end_atom_tag);
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

void MolModel::clear()
{
    QString desc = QString("Clear");
    auto redo = [this]() { clearFromCommand(); };
    doCommandWithMolUndo(redo, desc);
}

int MolModel::getTagForAtom(const RDKit::Atom* const atom)
{
    for (auto& [bookmark, atom_list] : *m_mol.getAtomBookmarks()) {
        // atom_list should be exactly one element long since we ensure
        // uniqueness for atom tags
        if (atom == atom_list.front()) {
            return bookmark;
        }
    }
    throw std::runtime_error("No atom tag found");
}

int MolModel::getTagForBond(const RDKit::Bond* const bond)
{
    for (auto& [bookmark, bond_list] : *m_mol.getBondBookmarks()) {
        // bond_list should be exactly one element long since we ensure
        // uniqueness for bond tags
        if (bond == bond_list.front()) {
            return bookmark;
        }
    }
    throw std::runtime_error("No bond tag found");
}

void MolModel::addAtomFromCommand(const int atom_tag,
                                  const std::string& element,
                                  const RDGeom::Point3D& coords)
{
    Q_ASSERT(m_in_command);
    // RDKit will take ownership of this atom after we call addAtom, so we don't
    // have to worry about deleting it
    auto* atom = new RDKit::Atom(element);
    // addAtom returns the atom index of the newly added atom (*not* the new
    // number of atoms)
    unsigned int atom_index = m_mol.addAtom(atom, /* updateLabel = */ false,
                                            /* takeOwnership = */ true);
    m_mol.setAtomBookmark(atom, atom_tag);
    m_mol.getConformer().setAtomPos(atom_index, coords);
    emit moleculeChanged();
}

void MolModel::removeAtomFromCommand(const int atom_tag)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* atom = m_mol.getUniqueAtomWithBookmark(atom_tag);
    m_mol.removeAtom(atom);
    emit moleculeChanged();
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
    m_mol.setBondBookmark(bond, bond_tag);
    emit moleculeChanged();
}

void MolModel::removeBondFromCommand(const int bond_tag,
                                     const int start_atom_tag,
                                     const int end_atom_tag)
{
    Q_ASSERT(m_in_command);
    RDKit::Atom* start_atom = m_mol.getUniqueAtomWithBookmark(start_atom_tag);
    RDKit::Atom* end_atom = m_mol.getUniqueAtomWithBookmark(end_atom_tag);
    m_mol.removeBond(start_atom->getIdx(), end_atom->getIdx());
    emit moleculeChanged();
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
        m_mol.setAtomBookmark(atom, cur_atom_tag);
        ++atom_index;
    }

    for (int cur_bond_tag : bond_tags) {
        RDKit::Bond* bond = m_mol.getBondWithIdx(bond_index);
        m_mol.setBondBookmark(bond, cur_bond_tag);
        ++bond_index;
    }
    emit moleculeChanged();
}

void MolModel::clearFromCommand()
{
    Q_ASSERT(m_in_command);
    m_mol.clear();
    emit moleculeChanged();
}

} // namespace sketcher
} // namespace schrodinger

#include "schrodinger/sketcher/molviewer/mol_model.moc"
