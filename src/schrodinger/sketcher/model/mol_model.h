#pragma once

#include <functional>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_set>
#include <utility>

#include <GraphMol/RWMol.h>
#include <GraphMol/QueryAtom.h>
#include <GraphMol/QueryBond.h>
#include <QUndoStack>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/abstract_undoable_model.h"
#include "schrodinger/sketcher/model/sketcher_model.h"

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

enum class MergeId {
    NO_MERGE = -1,
    ROTATE = 1,
    TRANSLATE,
};

/**
 * A model for making undoable changes to an RDKit Mol using a QUndoStack.
 * Note that all public methods in this class should fall into one of two
 * categories:
 *
 * Getters: These methods return a value but do not modify the instance
 * (i.e. these methods should be marked as const).
 *
 * Undoable commands: These methods undoably modify the instance by creating
 * a command and pushing it onto the undo stack (which immediately runs the
 * command) using AbstractUndoableModel::doCommand.  Note that these methods
 * should have a void return type, as callers should not be taking *any*
 * action in response to these methods being called.  Callers should instead
 * listen for signals, as this ensures that the caller will respond
 * identically when a command is initially done and when it's redone.
 *
 * In other words, don't modify addAtom to return the newly added atom. That
 * encourages the calling class to write code like
 *
 *   Atom* new_atom = mol_model.addAtom("C", coord);
 *   addGraphicsItemForAtom(new_atom);
 *
 * which is wrong!  addGraphicsItemForAtom will be called when the atom is
 * initially added, but it will *not* be called if the user undoes and
 * redoes the addAtom call.  Instead, the calling class should listen for
 * the moleculeChanged signal.
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
     * Undoably add a single atom.  This atom may optionally be bound to an
     * existing atom.
     *
     * @param element The element for the new atom
     * @param coords The coordinates for the new atom
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @param bond_type If bound_to is given, the type of bond that will be
     * added between the new atom and bound_to
     * @param bond_dir If bound_to is given, the stereochemistry of the bond
     * that will be added between the new atom and bound_to
     */
    void addAtom(
        const Element& element, const RDGeom::Point3D& coords,
        const RDKit::Atom* const bound_to_atom = nullptr,
        const RDKit::Bond::BondType& bond_type = RDKit::Bond::BondType::SINGLE,
        const RDKit::Bond::BondDir& bond_dir = RDKit::Bond::BondDir::NONE);

    /**
     * Undoably add a single query atom.  This atom may optionally be bound to
     * an existing atom.
     *
     * @param atom_query The query for the atom
     * @param coords The coordinates for the new atom
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @param bond_type If bound_to is given, the type of bond that will be
     * added between the new atom and bound_to
     * @param bond_dir If bound_to is given, the stereochemistry of the bond
     * that will be added between the new atom and bound_to
     */
    void addAtom(
        const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
        const RDGeom::Point3D& coords,
        const RDKit::Atom* const bound_to_atom = nullptr,
        const RDKit::Bond::BondType& bond_type = RDKit::Bond::BondType::SINGLE,
        const RDKit::Bond::BondDir& bond_dir = RDKit::Bond::BondDir::NONE);

    /**
     * Undoably add a single atom and bond it to an existing atom using a query
     * bond.
     *
     * @param element The element for the new atom
     * @param coords The coordinates for the new atom
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @param bond_query The query for the bond
     */
    void addAtom(
        const Element& element, const RDGeom::Point3D& coords,
        const RDKit::Atom* const bound_to_atom,
        const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query);

    /**
     * Undoably add a single query atom and bond it to an existing atom using a
     * query bond.
     *
     * @param atom_query The query for the atom
     * @param coords The coordinates for the new atom
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @param bond_query The query for the bond
     */
    void addAtom(
        const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
        const RDGeom::Point3D& coords, const RDKit::Atom* const bound_to_atom,
        const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query);

    /**
     * Undoably add an R-group atom.  This atom may optionally be bound to an
     * existing atom.
     *
     * @param r_group_num The R-group number for the new atom
     * @param coords The coordinates for the new atom
     * @param bound_to_atom If not nullptr, a single bond will be added between
     * this atom and the newly added atom
     */
    void addRGroup(const unsigned int r_group_num,
                   const RDGeom::Point3D& coords,
                   const RDKit::Atom* const bound_to_atom = nullptr);

    /**
     * Undoably add an attachment point bound to an existing atom.
     *
     * @param coords The coordinates for the new attachment point
     * @param bound_to_atom A single bond will be added between
     * this atom and the newly added atom.  Must not be nullptr.
     */
    void addAttachmentPoint(const RDGeom::Point3D& coords,
                            const RDKit::Atom* const bound_to_atom);

    /**
     * Undoably add a chain of atoms, where each atom is bound to the previous
     * and next atoms in the chain.
     *
     * @param element The element for the new atoms
     * @param coords The coordinates for the new atoms
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @param bond_type The type of bond to add
     * @param bond_dir The stereochemistry of the bond to add
     */
    void addAtomChain(
        const Element& element, const std::vector<RDGeom::Point3D>& coords,
        const RDKit::Atom* const bound_to_atom = nullptr,
        const RDKit::Bond::BondType& bond_type = RDKit::Bond::BondType::SINGLE,
        const RDKit::Bond::BondDir& bond_dir = RDKit::Bond::BondDir::NONE);

    /**
     * Undoably add a chain of query atoms, where each query atom is bound to
     * the previous and next query atoms in the chain.
     *
     * @param atom_query The query for the atoms
     * @param coords The coordinates for the new atoms
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @param bond_type The type of bond to add
     * @param bond_dir The stereochemistry of the bond to add
     */
    void addAtomChain(
        const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
        const std::vector<RDGeom::Point3D>& coords,
        const RDKit::Atom* const bound_to_atom = nullptr,
        const RDKit::Bond::BondType& bond_type = RDKit::Bond::BondType::SINGLE,
        const RDKit::Bond::BondDir& bond_dir = RDKit::Bond::BondDir::NONE);

    /**
     * Undoably add a chain of atoms, where each atom is bound to the previous
     * and next atoms in the chain using a query bond.
     *
     * @param element The element for the new atoms
     * @param coords The coordinates for the new atoms
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @param bond_query The query for the bonds
     */
    void addAtomChain(
        const Element& element, const std::vector<RDGeom::Point3D>& coords,
        const RDKit::Atom* const bound_to_atom,
        const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query);

    /**
     * Undoably add a chain of query atoms, where each query atom is bound to
     * the previous and next query atoms in the chain using a query bond.
     *
     * @param atom_query The query for the atoms
     * @param coords The coordinates for the new atoms
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @param bond_query The query for the bonds
     */
    void addAtomChain(
        const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query,
        const std::vector<RDGeom::Point3D>& coords,
        const RDKit::Atom* const bound_to_atom,
        const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query);

    /**
     * Undoably add a chain of R-group atoms, where each R-group atom is bound
     * to the previous and next atoms in the chain using a single bond.
     *
     * @param r_group_nums The R-group number for each new atom
     * @param coords The coordinates for the new atoms
     * @param bound_to_atom If not nullptr, a single bond will be added between
     * this atom and the first atom in the chain
     */
    void addRGroupChain(const std::vector<unsigned int> r_group_nums,
                        const std::vector<RDGeom::Point3D>& coords,
                        const RDKit::Atom* const bound_to_atom = nullptr);

    /**
     * Undoably add a bond between the specified atoms.
     *
     * @param start_atom The start atom for the bond
     * @param end_atom The end atom for the bond
     * @param bond_type The type of bond to add
     * @param bond_dir The stereochemistry of the bond to add
     */
    void
    addBond(const RDKit::Atom* const start_atom,
            const RDKit::Atom* const end_atom,
            const RDKit::Bond::BondType& bond_type,
            const RDKit::Bond::BondDir& bond_dir = RDKit::Bond::BondDir::NONE);

    /**
     * Undoably add a query bond between the specified atoms.
     *
     * @param start_atom The start atom for the bond
     * @param end_atom The end atom for the bond
     * @param bond_query The query for the bond
     */
    void addBond(
        const RDKit::Atom* const start_atom, const RDKit::Atom* const end_atom,
        const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query);

    /**
     * Undoably remove the given atoms and bonds.  Note that, even though this
     * method accepts pointers to const Atoms and Bonds, the Atoms and Bonds
     * will be destroyed as a result of this method and will no longer be valid.
     * (Without the const, it wouldn't be possible to pass in values returned
     * from getMol.)
     *
     * @note If either an attachment point dummy atom or the associated bond are
     * removed, then the other will automatically be removed as well.
     */
    void
    removeAtomsAndBonds(const std::unordered_set<const RDKit::Atom*>& atoms,
                        const std::unordered_set<const RDKit::Bond*>& bonds);

    /**
     * Unodable flip all atoms coordinates horizontally or vertically.
     */
    void flipAllHorizontal();
    void flipAllVertical();

    /**
     * rotate all atoms by the given angle (the rdkit coordinates rotate
     * clockwise, while the representation in the scene
     * rotates counter-clockwise,  since the y axis is inverted)
     * @param angle angle in degrees
     * @param pivot_point point to rotate around
     * @param atoms atoms to rotate, if empty rotate all atoms
     */
    void rotateByAngle(float angle, const RDGeom::Point3D& pivot_point,
                       const std::vector<const RDKit::Atom*>& atoms = {});

    /**
     * translate all atoms by the given vector
     * @param vector vector to translate by
     * @param atoms atoms to translate, if empty translate all atoms
     */
    void translateByVector(const RDGeom::Point3D& vector,
                           const std::vector<const RDKit::Atom*>& atoms = {});

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
     * Change the element of an existing atom.  This method can also mutate a
     * query atom into a non-query atom.
     */
    void mutateAtom(const RDKit::Atom* const atom, const Element& element);

    /**
     * Change the query of an existing atom.  This method can also mutate a
     * non-query atom into a query atom.
     */
    void mutateAtom(
        const RDKit::Atom* const atom,
        const std::shared_ptr<RDKit::QueryAtom::QUERYATOM_QUERY> atom_query);

    /**
     * Change the R-group of an existing atom.  This method can also mutate a
     * non-R-group atom into an R-group atom.
     */
    void mutateRGroup(const RDKit::Atom* const atom,
                      const unsigned int r_group_num);

    /**
     * Change the type of an existing bond.  This method can also mutate a
     * query bond into a non-query bond.
     *
     * @param bond The bond to mutate
     * @param bond_type The type of bond to mutate to
     * @param bond_dir The stereochemistry to mutate to
     */
    void mutateBond(
        const RDKit::Bond* const bond, const RDKit::Bond::BondType& bond_type,
        const RDKit::Bond::BondDir& bond_dir = RDKit::Bond::BondDir::NONE);

    /**
     * Change the query of an existing bond.  This method can also mutate a
     * non-query bond into a query bond.
     *
     * @param bond The bond to mutate
     * @param bond_query The query for the bond
     */
    void mutateBond(
        const RDKit::Bond* const bond,
        const std::shared_ptr<RDKit::QueryBond::QUERYBOND_QUERY> bond_query);

    /**
     * Reverse an existing bond (i.e. swap the start and end atoms).  This is
     * only meaningful for dative bonds or bonds with a BondDir of BEGINWEDGE or
     * BEGINDASH.
     */
    void flipBond(const RDKit::Bond* const bond);

    /**
     * Fully generate coordinates for the current molecule
     */
    void regenerateCoordinates();

    /**
     * Undoably clear the molecule.
     */
    void clear();

    /**
     * Undoably select or deselect the specified atoms and bonds.
     *
     * @param atoms The atoms to select or deselect
     * @param bonds The bonds to select or deselect
     * @param select_mode Whether to select, deselect, toggle selection, or
     * select-only (i.e. clear the selection and then select)
     *
     * @note If either an attachment point dummy atom or the associated bond are
     * selected (or deselected, etc.), then the other will automatically be
     * selected (or deselected, etc.) as well.
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
     * Signal emitted when the molecule is changed (e.g. by adding or removing
     * atoms or bonds), excluding selection changes, which will trigger
     * selectionChanged instead and coordinates changes which will trigger
     * coordinatesChanged.
     */
    void moleculeChanged();

    /**
     * Signal emitted when the molecule's coordinates are changed
     */
    void coordinatesChanged();

    /**
     * Signal emitted when selection is changed.  Note that atoms and bonds will
     * automatically be deselected when they are removed and reselected if the
     * removal is undone, and this signal *will* be emitted in those scenario.
     */
    void selectionChanged();

  protected:
    RDKit::RWMol m_mol = RDKit::RWMol();

    int m_next_atom_tag = 0;
    int m_next_bond_tag = 0;
    std::unordered_set<int> m_selected_atom_tags;
    std::unordered_set<int> m_selected_bond_tags;

    /**
     * Transform all coordinates using the given function.
     * @param desc The description to use for the redo/undo command.
     * @param function The function to use to transform the coordinates.
     * @param merge_id The merge id to use for the redo/undo command. If this is
     * different from -1, the command will be merged with the previous command
     * if they share the same merge id.
     * @param atoms The atoms to transform. If this is empty, all atoms will be
     * transformed
     */
    void transformCoordinatesWithFunction(
        const QString& desc, std::function<void(RDGeom::Point3D&)> function,
        MergeId merge_id = MergeId::NO_MERGE,
        std::vector<const RDKit::Atom*> atoms =
            std::vector<const RDKit::Atom*>());

    /**
     * Set the atom tag for the specified atom
     */
    void setTagForAtom(RDKit::Atom* const atom, const int atom_tag);

    /**
     * Find the atom tag for the specified atom.
     *
     * Note that this method is const, but not marked as such because RDKit
     * does not implement a const version of getAtomBookmarks.  (And because
     * this method isn't public, so it's not worth using const_cast unless we
     * need this to be const.)
     */
    int getTagForAtom(const RDKit::Atom* const atom);

    /**
     * Set the bond tag for the specified bond
     */
    void setTagForBond(RDKit::Bond* const bond, const int bond_tag);

    /**
     * Find the bond tag for the specified bond.
     *
     * Note that this method is const, but not marked as such because RDKit
     * does not implement a const version of getBondBookmarks.  (And because
     * this method isn't public, so it's not worth using const_cast unless we
     * need this to be const.)
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
     * Determine all atom and bond tags needed for adding an atom chain
     * @param coords The list of coordinates for the atom chain being added
     * @param bound_to_atom If not nullptr, a bond will be added between this
     * atom and the first atom in the chain
     * @return A tuple of
     *   - The atom tags
     *   - The bond tags
     *   - The tag for bound_to_atom (or -1 if bound_to_atom was nullptr)
     */
    std::tuple<std::vector<int>, std::vector<int>, int>
    getAtomAndBondTagsForAddingAtomChain(
        const std::vector<RDGeom::Point3D>& coords,
        const RDKit::Atom* const bound_to_atom);

    /**
     * Get the atomic number and element name for the specified element.
     */
    std::pair<unsigned int, QString> getAddElementInfo(const Element& element);

    /**
     * Undoably select or deselect the specified atoms and bonds.
     *
     * @param atom_tags Tags for the atoms to select or deselect
     * @param bond_tags Tags for the bonds to select or deselect
     * @param select_mode Whether to select, deselect, toggle selection, or
     * select-only (i.e. clear the selection and then select)
     *
     * @note This method does *not* expand attachment points (i.e. make sure
     * that both the dummy atom and the bond are included).  That is instead
     * handled in the select method.
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
     * Undoably select or deselect the specified atoms and bonds.  If
     * selecting, all given tags must be currently unselected.  If deselecting,
     * all given tags must be currently selected.
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
     * Update all molecule metadata and notify the Scene of the changes.  This
     * method should be called exactly once after changes to m_mol.
     */
    void finalizeMoleculeChange(bool selection_changed = false);

    /**
     * Update any RDKit metadata that is required to render the molecule in the
     * Scene.  This method must be called anytime atoms or bonds are changed.
     */
    void updateMoleculeMetadata();

    /**
     * Generate multiple atom tags or bond tags
     * @param[in] count The number of tags to generate
     * @param[in,out] tag_counter Either m_next_atom_tag or m_next_bond_tag,
     * depending on whether atom tags or bond tags are desired.
     * @return The generated tags
     */
    std::vector<int> getNextNTags(const size_t count, int& tag_counter) const;

    /**
     * Make sure that the returned sets include both an attachment point dummy
     * atom *and* the associated bond if either the atom *or* the bond are
     * included in the input sets.
     */
    std::pair<std::unordered_set<const RDKit::Atom*>,
              std::unordered_set<const RDKit::Bond*>>
    ensureCompleteAttachmentPoints(
        const std::unordered_set<const RDKit::Atom*>& atoms,
        const std::unordered_set<const RDKit::Bond*>& bonds);

    /**
     * Add a chain of atoms, where each atom is bound to the previous and next
     * atoms in the chain.  This method must only be called from an undo
     * command.
     *
     * @param atom_tags The atom tags for each newly created atom
     * @param bond_tags The bond tags for each newly created bond
     * @param create_atom A function that will return a new atom to add.  Note
     * that this method will add a *copy* of the returned atom, as the
     * shared_ptr is responsible for deleting the returned atom.
     * @param coords The coordinates for the new atoms
     * @param create_bond A function that will return a new bond to add.  Note
     * that this method will add a *copy* of the returned bond, as the
     * shared_ptr is responsible for deleting the returned bond.
     * @param bound_to_atom_tag If given, a bond will be added between the
     * existing atom with this tag and the first new atom in the chain
     */
    void addAtomChainFromCommand(
        const std::vector<int>& atom_tags, const std::vector<int>& bond_tags,
        const std::function<std::shared_ptr<RDKit::Atom>()> create_atom,
        const std::vector<RDGeom::Point3D>& coords,
        const std::function<std::shared_ptr<RDKit::Bond>()> create_bond,
        const int bound_to_atom_tag);

    /**
     * Remove an atom from the molecule.  This method must only be called from
     * an undo command.  Note that this method does not update ring information
     * or emit signals, as these are intended to be done from
     * removeAtomsAndBondsFromCommand instead.
     * @param atom_tag The tag of the atom to delete
     * @return whether selection was changed by this action
     */
    bool removeAtomFromCommand(const int atom_tag);

    /**
     * Add a bond to the molecule.  This method must only be called from an undo
     * command.
     * @param bond_tag The tag to use for the newly added bond
     * @param start_atom_tag The tag of the bond's start atom
     * @param end_atom_tag The tag of the bond's end atom
     * @param create_bond A function that will return the new bond to add.  Note
     * that this method will add a *copy* of the returned bond, as the
     * shared_ptr is responsible for deleting the returned bond.
     */
    void addBondFromCommand(
        const int bond_tag, const int start_atom_tag, const int end_atom_tag,
        const std::function<std::shared_ptr<RDKit::Bond>()> create_bond);

    /**
     * Remove a bond from the molecule.  This method must only be called from
     * an undo command.  Note that this method does not update ring information
     * or emit signals (as these are intended to be done from
     * removeAtomsAndBondsFromCommand instead).
     * @return whether selection was changed by this action
     */
    bool removeBondFromCommand(const int bond_tag, const int start_atom_tag,
                               const int end_atom_tag);

    /**
     * set new coordinates for a set of atoms.  This method must only be called
     * from an undo command.
     */
    void setCoordinates(const std::vector<int>& atom_tags,
                        const std::vector<RDGeom::Point3D>& coords);

    /**
     * Remove the specified atoms and bonds from the molecule.  This method
     * must only be called from an undo command.
     * @param atom_tags The atom tags to delete
     * @param bond_tags_with_atoms A list of (bond tag, bond's start atom tag,
     * bond's end atom tag) for the bonds to delete
     */
    void removeAtomsAndBondsFromCommand(
        const std::vector<int>& atom_tags,
        const std::vector<std::tuple<int, int, int>>& bond_tags_with_atoms);

    /**
     * Add all atoms and bonds from the given molecule into this molecule. This
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
     * Change the element of an existing atom, or change a query atom to a
     * non-query atom.  This method must only be called from an undo command.
     * @param atom_tag The tag of the atom to mutate
     * @param create_atom A function that will return a new atom to mutate to.
     * Note that this method will add a *copy* of the returned atom, as the
     * shared_ptr is responsible for deleting the returned atom.
     */
    void mutateAtomFromCommand(
        const int atom_tag,
        const std::function<std::shared_ptr<RDKit::Atom>()> create_atom);

    /**
     * Change the type of an existing bond, or change a query bond to a
     * non-query bond.  This method must only be called from an undo command.
     * @param bond_tag The tag of the bond to mutate
     * @param create_bond A function that will return a new bond to mutate to.
     * Note that this method will add a *copy* of the returned bond, as the
     * shared_ptr is responsible for deleting the returned bond.
     */
    void mutateBondFromCommand(
        const int bond_tag,
        const std::function<std::shared_ptr<RDKit::Bond>()> create_bond);

    /**
     * Reverse an existing bond (i.e. swap the start and end atoms).  This
     * method must only be called from an undo command.
     */
    void flipBondFromCommand(const int bond_tag);

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
     * Clear all selected atoms and bonds.  This method must only be called
     * from an undo command.
     */
    void clearSelectionFromCommand();
};

} // namespace sketcher
} // namespace schrodinger
