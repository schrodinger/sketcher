/**
 * Copyright Schrodinger, LLC. All rights reserved.
 */
#pragma once

#include <functional>
#include <string>

#include <GraphMol/RWMol.h>

#include <QUndoStack>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_undoable_model.h"

class QPointF;
class QString;
class QUndoStack;

namespace RDKit
{
class ROMol;
class Atom;
class Bond;
} // namespace RDKit

namespace schrodinger
{
namespace sketcher
{

/**
 * A model for making undoable changes to an RDKit Mol using a QUndoStack.  Note
 * that all public methods in this class should fall into one of two categories:
 *
 * Getters: These methods return a value but do not modify the instance (i.e.
 * these methods should be marked as const).
 *
 * Undoable commands: These methods undoably modify the instance by creating a
 * command and pushing it onto the undo stack (which immediately runs the
 * command) using AbstractUndoableModel::doCommand.  Note that these methods
 * should have a void return type, as callers should not be taking *any* action
 * in response to these methods being called.  Callers should instead listen for
 * signals, as this ensures that the caller will respond identically when a
 * command is initially done and when it's redone.
 *
 * In other words, don't modify addAtom to return the newly added atom.  That
 * encourages the calling class to write code like
 *
 *   Atom* new_atom = mol_model.addAtom("C", coord);
 *   addGraphicsItemForAtom(new_atom);
 *
 * which is wrong!  addGraphicsItemForAtom will be called when the atom is
 * initially added, but it will *not* be called if the user undoes and redoes
 * the addAtom call.  Instead, the calling class should listen for the
 * moleculeChanged signal.
 */
class SKETCHER_API MolModel : public AbstractUndoableModel
{
    Q_OBJECT
  public:
    MolModel(QUndoStack* const undo_stack = nullptr);

    /******************************** GETTERS *******************************/

    /**
     * @return the RDKit molecule
     */
    const RDKit::ROMol* getMol() const;

    /*************************** UNDOABLE COMMANDS **************************/

    /**
     * Undoably add the specified atom.
     *
     * @param element The element for the new atom
     * @param coords The coordinates for the new atom
     */
    void addAtom(const std::string& element, const RDGeom::Point3D& coords);
    // TODO: accept element constant here?
    // TODO: add support for query atoms, R groups

    /**
     * Undoably remove the given atom.  Note that even though this method
     * accepts a pointer to a const Atom, the Atom will be destroyed as a result
     * of this method and will no longer be valid.  (Without the const, it
     * wouldn't be possible to pass in a value returned from getMol.)
     */
    void removeAtom(const RDKit::Atom* const atom);

    /**
     * Undoably add a bond between the specified atoms.
     */
    void addBond(const RDKit::Atom* const start_atom,
                 const RDKit::Atom* const end_atom,
                 const RDKit::Bond::BondType& bond_type);
    // TODO: add support for query bonds

    /**
     * Undoably remove the given bond.  Note that, even though this method
     * accepts a pointer to a const Bond, the Bond will be destroyed as a result
     * of this method and will no longer be valid.  (Without the const, it
     * wouldn't be possible to pass in a value returned from getMol.)
     */
    void removeBond(const RDKit::Bond* const bond);

    /**
     * Undoably add all atoms and bonds from the given molecule into this
     * molecule.
     *
     * @param mol The molecule to add
     * @param description The description to use for the undo command.
     */
    void addMol(const RDKit::ROMol& mol,
                const QString& description = "Import molecule");

    /**
     * Undoably clear the molecule.
     */
    void clear();

    // TODO:
    // void modifyAtom(/*arguments here*/);
    // void modifyBond(/*arguments here*/);
    // void clearSelection();
    // void select(/*arguments here*/);

  signals:

    /**
     * Signal emitted when the molecule is changed, excluding selection changes,
     * which will trigger selectionChanged instead.
     */
    void moleculeChanged();

    // TODO:
    void selectionChanged(/*params here*/);
    void selectionCleared();

  protected:
    RDKit::RWMol m_mol = RDKit::RWMol();
    int m_next_atom_tag = 0;
    int m_next_bond_tag = 0;

    /**
     * Find the atom tag for the specified atom.
     *
     * Note that this method is const, but not marked as such because RDKit does
     * not implement a const version of getAtomBookmarks.  (And because this
     * method isn't public, so it's not worth using const_cast unless we need
     * this to be const.)
     *
     * @throw std::runtime_error if the atom is not part of this molecule
     */
    int getTagForAtom(const RDKit::Atom* const atom);

    /**
     * Find the bond tag for the specified bond.
     *
     * Note that this method is const, but not marked as such because RDKit does
     * not implement a const version of getBondBookmarks.  (And because this
     * method isn't public, so it's not worth using const_cast unless we need
     * this to be const.)
     *
     * @throw std::runtime_error if the bond is not part of this molecule
     */
    int getTagForBond(const RDKit::Bond* const bond);

    /**
     * Create an undo command and add it to the undo stack, which automatically
     * executes the command.  The undo operation for this command will be
     * carried out by restoring the entire RDKit molecule.  This preserves
     * atom/bond indices and all atom/bond properties, but it means that an
     * undo_stack.undo() call will invalidate all RDKit Atom and Bond objects.
     *
     * @param redo The function to call on redo (or for the initial do)
     * @param description A description of the command
     */
    void doCommandWithMolUndo(const std::function<void()> redo,
                              const QString& description);

    /**
     * Add an atom to the molecule.  This method must only be called from an
     * undo command.
     */
    void addAtomFromCommand(const int atom_tag, const std::string& element,
                            const RDGeom::Point3D& coords);

    /**
     * Remove an atom from the molecule.  This method must only be called from
     * an undo command.
     */
    void removeAtomFromCommand(const int atom_tag);

    /**
     * Add a bond to the molecule.  This method must only be called from an undo
     * command.
     */
    void addBondFromCommand(const int bond_tag, const int start_atom_tag,
                            const int end_atom_tag,
                            const RDKit::Bond::BondType& bond_type);

    /**
     * Remove a bond from the molecule.  This method must only be called from an
     * undo command.
     */
    void removeBondFromCommand(const int bond_tag, const int start_atom_tag,
                               const int end_atom_tag);

    /**
     * Add all atoms and bonds from the given molecule into this molecule.  This
     * method must only be called from an undo command.
     *
     * @param mol The molecule to add
     * @param atom_tags A list of atom tags to use for each new atom.  The
     * length of this list must be exactly mol.getNumAtoms().
     * @param bond_tags A list of bond tags to use for each new bond.  The
     * length of this list must be exactly mol.getNumBonds().
     */
    void addMolFromCommand(const RDKit::ROMol& mol,
                           const std::vector<int>& atom_tags,
                           const std::vector<int>& bond_tags);

    /**
     * Clear the molecule.  This method must only be called from an undo
     * command.
     */
    void clearFromCommand();
};

} // namespace sketcher
} // namespace schrodinger
