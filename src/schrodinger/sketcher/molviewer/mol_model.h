#pragma once

#include <functional>
#include <string>
#include <unordered_set>
#include <utility>

#include <GraphMol/RWMol.h>
#include <QUndoStack>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/molviewer/abstract_undoable_model.h"

class QObject;
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

namespace rdkit_extensions
{
enum class Format;
}

namespace sketcher
{

enum class SelectMode {
    SELECT,      // Shift + click
    DESELECT,    //
    TOGGLE,      // Ctrl + click
    SELECT_ONLY, // normal click (no Shift or Ctrl)
};

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
    MolModel(QUndoStack* const undo_stack = nullptr, QObject* parent = nullptr);

    /******************************** GETTERS *******************************/

    /**
     * @return the RDKit molecule
     */
    const RDKit::ROMol* getMol() const;

    /**
     * @param format format to convert to
     * @return serialized representation of the sketcher contents
     */
    std::string getMolText(rdkit_extensions::Format format);

    /**
     * @return A set of all currently selected atoms
     */
    std::unordered_set<const RDKit::Atom*> getSelectedAtoms() const;

    /**
     * @return A set of all currently selected bonds
     */
    std::unordered_set<const RDKit::Bond*> getSelectedBonds() const;

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
     * Create a molecule from a text string and load that into the scene.
     * Atomic coordinates will be automatically generated using coordgen.
     *
     * @param text input data to load
     * @param format format to parse
     */
    void addMolFromText(const std::string& text,
                        rdkit_extensions::Format format);

    /**
     * Undoably clear the molecule.
     */
    void clear();

    // TODO:
    // void modifyAtom(/*arguments here*/);
    // void modifyBond(/*arguments here*/);

    /**
     * Undoably select or deselect the specified atoms and bonds.
     *
     * @param atoms The atoms to select or deselect
     * @param bonds The bonds to select or deselect
     * @param select_mode Whether to select, deselect, toggle selection, or
     * select-only (i.e. clear the selection and then select)
     */
    void select(const std::unordered_set<const RDKit::Atom*>& atoms,
                const std::unordered_set<const RDKit::Bond*>& bonds,
                const SelectMode select_mode);

    /**
     * Undoably clear all selected atoms and bonds.
     */
    void clearSelection();

    /**
     * Select all atoms and bonds.
     */
    void selectAll();

    /**
     * Select all unselected atoms and bonds and deselect all selected atoms
     * and bonds.
     */
    void invertSelection();

  signals:

    /**
     * Signal emitted when the molecule is changed, excluding selection changes,
     * which will trigger selectionChanged instead.
     */
    void moleculeChanged();

    // Signal emitted when selection is changed.  Note that atoms and bonds will
    // automatically be deselected when they are removed and reselected if the
    // removal is undone, and this signal *will* be emitted in those scenario.
    void selectionChanged();

  protected:
    RDKit::RWMol m_mol = RDKit::RWMol();
    int m_next_atom_tag = 0;
    int m_next_bond_tag = 0;
    std::unordered_set<int> m_selected_atom_tags;
    std::unordered_set<int> m_selected_bond_tags;

    /**
     * Set the atom tag for the specified atom
     */
    void setTagForAtom(RDKit::Atom* const atom, const int atom_tag);

    /**
     * Find the atom tag for the specified atom.
     *
     * Note that this method is const, but not marked as such because RDKit does
     * not implement a const version of getAtomBookmarks.  (And because this
     * method isn't public, so it's not worth using const_cast unless we need
     * this to be const.)
     */
    int getTagForAtom(const RDKit::Atom* const atom);

    /**
     * Set the bond tag for the specified bond
     */
    void setTagForBond(RDKit::Bond* const bond, const int bond_tag);

    /**
     * Find the bond tag for the specified bond.
     *
     * Note that this method is const, but not marked as such because RDKit does
     * not implement a const version of getBondBookmarks.  (And because this
     * method isn't public, so it's not worth using const_cast unless we need
     * this to be const.)
     */
    int getTagForBond(const RDKit::Bond* const bond);

    /**
     * Return the atom identified by the given atom tag.  Note that an atom tag
     * uniquely identifies a single atom in the molecule.  No two atoms in the
     * same molecule will ever have identical atom tags, even if the two atoms
     * don't exist at the same time.
     *
     * @param atom_tag An atom_tag for an atom that is currently in the
     * molecule.
     *
     * @throw if no atom is found with the specified atom_tag
     */
    const RDKit::Atom* getAtomFromTag(int atom_tag) const;

    /**
     * Return the bond identified by the given bond tag.  Note that a bond tag
     * uniquely identifies a single bond in the molecule.  No two bonds in the
     * same molecule will ever have identical bond tags, even if the two bonds
     * don't exist at the same time.
     *
     * @param bond_tag A bond tag for a bond that is currently in the molecule.
     *
     * @throw if no bond is found with the specified bond_tag
     */
    const RDKit::Bond* getBondFromTag(int bond_tag) const;

    /**
     * Undoably select or deselect the specified atoms and bonds.
     *
     * @param atom_tags Tags for the atoms to select or deselect
     * @param bond_tags Tags for the bonds to select or deselect
     * @param select_mode Whether to select, deselect, toggle selection, or
     * select-only (i.e. clear the selection and then select)
     */
    void selectTags(const std::unordered_set<int>& atom_tags,
                    const std::unordered_set<int>& bond_tags,
                    const SelectMode select_mode);

    /**
     * Divide the given set of tags into two sets: one containing selected tags
     * and one containing unselected tags
     * @param tags_to_divide The atom or bond tags to divide
     * @param selected_tags The atom or bond tags that are currently selected
     * @return the selected and unselected sets, in that order
     */
    std::pair<std::unordered_set<int>, std::unordered_set<int>>
    divideBySelected(const std::unordered_set<int>& tags_to_divide,
                     const std::unordered_set<int>& selected_tags);

    /**
     * @return The atoms and bond tags for all atoms and bonds that are not
     * currently selected
     */
    std::pair<std::unordered_set<int>, std::unordered_set<int>>
    getAllUnselectedTags();

    /**
     * Undoably select or deselect the specified atoms and bonds.  If selecting,
     * all given tags must be currently unselected.  If deselecting, all given
     * tags must be currently selected.
     *
     * @param filtered_atom_tags Tags for the atoms to select or deselect.
     * @param filtered_bond_tags Tags for the bonds to select or deselect.
     * @param to_select Whether to select or deselect the specified atoms and
     * bonds.
     * @param description The description for the undo command
     */
    void doSelectionCommand(const std::unordered_set<int>& filtered_atom_tags,
                            const std::unordered_set<int>& filtered_bond_tags,
                            const bool to_select, const QString& description);

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

    /**
     * Select or deselect the specified atoms and bonds.  This method must only
     * be called from an undo command.
     *
     * @param atom_tags The atom tags to select or deselect
     * @param bond_tags The bond tags to select or deselect
     * @param selected Whether to select or deselect the specified atoms and
     * bonds
     */
    void setSelectedFromCommand(const std::unordered_set<int>& atom_tags,
                                const std::unordered_set<int>& bond_tags,
                                const bool selected);

    /**
     * Clear all selected atoms and bonds.  This method must only be called from
     * an undo command.
     */
    void clearSelectionFromCommand();
};

} // namespace sketcher
} // namespace schrodinger
