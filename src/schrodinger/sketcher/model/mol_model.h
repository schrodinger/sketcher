#pragma once

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <rdkit/GraphMol/ChemReactions/Reaction.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/QueryBond.h>
#include <QUndoStack>

#include "schrodinger/sketcher/definitions.h"
#include "schrodinger/sketcher/model/abstract_undoable_model.h"
#include "schrodinger/sketcher/model/non_molecular_object.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/sketcher/rdkit/fragment.h"

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

enum class ExplicitHActions {
    REMOVE,
    ADD,
    TOGGLE // add explicit Hs if any are missing, otherwise remove them all
};

/**
 * Values used as parameters to doCommandUsingSnapshots that indicate what types
 * of data will be changed by the command.
 */
typedef uint8_t WhatChangedType;
namespace WhatChanged
{
enum : WhatChangedType { // clang-format off
    NOTHING      = 0,
    MOLECULE     = 1 << 0,
    NON_MOL_OBJS = 1 << 1,
    ALL          = MOLECULE | NON_MOL_OBJS, // clang-format on
};
}

struct HighlightingInfo {
    HighlightingInfo(const std::unordered_set<int>& atom_tags,
                     const std::unordered_set<int>& bond_tags,
                     const QColor& color) :
        atom_tags(atom_tags),
        bond_tags(bond_tags),
        color(color){};
    std::unordered_set<int> atom_tags;
    std::unordered_set<int> bond_tags;
    QColor color;
};

/**
 * A copy of a MolModel state
 */
struct MolModelSnapshot {
    RDKit::RWMol m_mol;
    std::vector<NonMolecularObject> m_pluses;
    std::optional<NonMolecularObject> m_arrow;
    std::unordered_set<int> m_selected_atom_tags;
    std::unordered_set<int> m_selected_bond_tags;
    std::unordered_set<int> m_selected_non_molecular_tags;
    std::vector<HighlightingInfo> m_highlighting_info;

    MolModelSnapshot(const RDKit::RWMol& mol,
                     const std::vector<NonMolecularObject>& pluses,
                     const std::optional<NonMolecularObject>& arrow,
                     const std::unordered_set<int>& selected_atom_tags,
                     const std::unordered_set<int>& selected_bond_tags,
                     const std::unordered_set<int>& selected_non_molecular_tags,
                     const std::vector<HighlightingInfo>& m_highlighting_info);

    /**
     * @return Whether the selection in the other snapshot is identical to the
     * selection in this snapshot.
     */
    bool isSelectionIdentical(const MolModelSnapshot& other);
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
     * @return a pointer to the actual RDKit molecule (i.e. not a copy)
     */
    const RDKit::ROMol* getMol() const;

    /**
     * @return the copy of the RDKit molecule with all internal MolModel
     * properties removed
     */
    boost::shared_ptr<RDKit::ROMol> getMolForExport() const;

    /**
     * @return the current selection as a new RDKit molecule with all internal
     * MolModel properties removed
     */
    boost::shared_ptr<RDKit::ROMol> getSelectedMolForExport();

    /**
     * @return The RDKit reaction contained in this model.  The reaction
     * molecules will have all internal MolModel properties removed.
     * @throw std::runtime_error if no reaction arrow is present
     */
    boost::shared_ptr<RDKit::ChemicalReaction> getReactionForExport() const;

    /**
     * @return whether the model is completely empty, i.e. no atoms, bonds, or
     * non-molecular objects
     */
    bool isEmpty() const;

    /**
     * @return whether anything (atoms, bonds, or non-molecular objects) is
     * selected
     */
    bool hasSelection() const;

    /**
     * @return halo highlighting information
     */
    std::vector<std::tuple<std::unordered_set<const RDKit::Atom*>,
                           std::unordered_set<const RDKit::Bond*>, QColor>>
    getHaloHighlighting() const;

    /**
     * clear all halo highlighting
     */
    void clearHaloHighlighting();

    /**
     * add halo highlighting
     * @param atoms atoms to highlight
     * @param bonds bonds to highlight
     * @param color color to highlight with
     */
    void addHaloHighlighting(const std::unordered_set<int>& atoms,
                             const std::unordered_set<int>& bonds,
                             const QColor& color);

    /**
     * @return A set of all currently selected atoms
     */
    std::unordered_set<const RDKit::Atom*> getSelectedAtoms() const;

    /**
     * @return A set of all currently selected bonds
     */
    std::unordered_set<const RDKit::Bond*> getSelectedBonds() const;

    /**
     * @return A set of all currently selected non-molecular objects
     */
    std::unordered_set<const NonMolecularObject*>
    getSelectedNonMolecularObjects() const;

    /**
     * @return whether the model contains a reaction arrow
     */
    bool hasReactionArrow() const;

    /**
     * @return the reaction arrow.  Will return nullptr if no reaction arrow is
     * present.
     */
    const NonMolecularObject* getReactionArrow() const;

    /**
     * @return all non-molecular objects (pluses and the reaction arrow)
     */
    std::unordered_set<const NonMolecularObject*>
    getNonMolecularObjects() const;

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
     * Undoably add a non-molecular object (a plus sign or a reaction arrow)
     * @param type The type of non-molecular object to add
     * @param coords The coordinates for the new non-molecular object
     */
    void addNonMolecularObject(const NonMolecularType& type,
                               const RDGeom::Point3D& coords);

    /**
     * Undoably remove the given atoms, bonds, sgroups and non-molecular objects
     * (pluses and/or the reaction arrow).  Note that, even though this method
     * accepts pointers to const objects, the passed in objects will be
     * destroyed as a result of this method and will no longer be valid.
     * (Without the const, it wouldn't be possible to pass in values returned
     * from getMol.)
     *
     * @note If either an attachment point dummy atom or the associated bond are
     * removed, then the other will automatically be removed as well.
     */
    void remove(const std::unordered_set<const RDKit::Atom*>& atoms,
                const std::unordered_set<const RDKit::Bond*>& bonds,
                const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
                const std::unordered_set<const NonMolecularObject*>&
                    non_molecular_objects);

    /**
     * Undoably flip the given bond, mirroring the coordinates of the
     * substituent identified by it and leaving the core of the molecule
     * unchanged. Substituent and core are defined as all the atoms and
     * bonds connected to one of the bond's atoms, excluding the bond itself.
     * The substituent is the set with the smaller number of atoms.
     */
    void flipSubstituent(const RDKit::Bond* const bond);

    /**
     * Undoably flip all atoms coordinates around the given axis.
     * @param p1 first point on the axis
     * @param p2 second point on the axis
     * @param atoms atoms to flip
     */
    void flipAroundSegment(const RDGeom::Point3D& p1, const RDGeom::Point3D& p2,
                           const std::unordered_set<const RDKit::Atom*>& atoms);

    /**
     * Undoably flip all atoms coordinates horizontally or vertically.
     */
    void flipAllHorizontal();
    void flipAllVertical();

    /**
     * rotate all atoms by the given angle (the rdkit coordinates rotate
     * clockwise, while the representation in the scene
     * rotates counter-clockwise,  since the y axis is inverted)
     * @param angle angle in degrees
     * @param pivot_point point to rotate around
     * @param atoms atoms to rotate.  If both atoms and non_molecular_objects
     * are empty, then all objects will be rotated.
     * @param non_molecular_objects non-molecular objects to rotate.  If both
     * atoms and non_molecular_objects are empty, then all objects will be
     * rotated.
     */
    void rotateByAngle(float angle, const RDGeom::Point3D& pivot_point,
                       const std::unordered_set<const RDKit::Atom*>& atoms = {},
                       const std::unordered_set<const NonMolecularObject*>&
                           non_molecular_objects = {});

    /**
     * translate all atoms by the given vector
     * @param vector vector to translate by
     * @param atoms atoms to translate.  If both atoms and non_molecular_objects
     * are empty, then all objects will be translated.
     * @param non_molecular_objects non-molecular objects to translate.  If both
     * atoms and non_molecular_objects are empty, then all objects will be
     * translated.
     */
    void
    translateByVector(const RDGeom::Point3D& vector,
                      const std::unordered_set<const RDKit::Atom*>& atoms = {},
                      const std::unordered_set<const NonMolecularObject*>&
                          non_molecular_objects = {});

    /**
     * Combine pairs of atoms into single atoms with all bonds from both
     * original atoms.  This is typically used to merge overlapping atoms at the
     * end of a rotation or translation.
     *
     * @param atoms_to_merge A vector of {atom to keep, stationary atom} pairs.
     * The resulting merged atom has the coordinates of the stationary atom, but
     * all other properties (e.g. element, charge) are taken from the atom to
     * keep.
     */
    void mergeAtoms(
        const std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>>&
            atoms_to_merge);

    /**
     * Undoably add all atoms and bonds from the given molecule into this
     * molecule. The molecule will be centered at the origin.
     *
     * @param mol The molecule to add
     * @param description The description to use for the undo command.
     * @param reposition_mol Whether to reposition the molecule.  The there is
     * an existing molecular structure in this model, the new molecule will be
     * positioned to the right of it.  Otherwise, the new molecule will be
     * centered.
     */
    void addMol(RDKit::RWMol mol,
                const QString& description = "Import molecule",
                const bool reposition_mol = true);

    /**
     * Undoably add the given reaction to the model.  All reactants and products
     * from the reaction will be added, and plus signs and arrows will be added
     * between the molecules.  Note that reaction agents will *not* be added.
     *
     * @throw std::runtime_error if the model already contains a reaction arrow
     */
    void addReaction(const RDKit::ChemicalReaction& reaction);

    /**
     * Undoably add a fragment and (optionally) bond it to the existing
     * structure.  The fragment must contain exactly one attachment point and
     * must have a conformer set.
     *
     * @param fragment_to_add The fragment to add
     * @param core_start_atom An existing atom with the same coordinates as the
     * fragment's attachment point parent atom.  If given, existing atoms may be
     * replaced by fragment atoms with the same coordinates, following the rules
     * given below.
     *
     * To determine which core atoms will be replaced, first we must determine
     * corresponding atom.  In order to "correspond," a fragment atom and an
     * existing atom must have the same coordinates.  Additionally, all
     * corresponding atoms must be connected.  In other words, frag_start_atom
     * always corresponds to core_start_atom.  Neighbors of frag_start_atom may
     * correspond to neighbors of core_start_atom, assuming they have the same
     * coordinates.  If a pair of those atoms are corresponding, then any
     * neighbors of *those* atoms may correspond if they have the same
     * coordinates, etc.
     *
     * If an existing atom has a corresponding fragment atom, then the core atom
     * will be replaced if the fragment atom is a query atom or a heteroatom, or
     * if the existing atom is neither a query atom nor a heteroatom.  All
     * fragment bonds will be added to the molecule, and any existing bonds
     * between corresponding atoms will be replaced by the corresponding
     * fragment bond.
     */
    void addFragment(const RDKit::ROMol& fragment_to_add,
                     const RDKit::Atom* const core_start_atom = nullptr);

    /**
     * Set explicit hydrogens on or off for the specified atoms
     * @param atoms The atoms to consider. If empty, all atoms will be used.
     * @param action wether to add or remove explicit Hs
     */
    void updateExplicitHs(const ExplicitHActions action,
                          std::unordered_set<const RDKit::Atom*> atoms = {});

    /**
     * Change the charge of an existing atom.
     * @param atom The atom to mutate
     * @param charge The charge to be set
     */
    void setAtomCharge(const RDKit::Atom* const atom, int charge);

    /**
     * Add mapping info to a set of existing atoms.
     * @param atoms The atoms to map
     * @param mapping The mapping number to be set. Setting 0 will remove the
     * mapping
     */
    void setAtomMapping(const std::unordered_set<const RDKit::Atom*>& atoms,
                        int mapping);

    /**
     * Reverse an existing bond (i.e. swap the start and end atoms).  This
     * is only meaningful for dative bonds or bonds with a BondDir of
     * BEGINWEDGE or BEGINDASH.
     */
    void flipBond(const RDKit::Bond* const bond);

    /**
     * Fully generate coordinates for the current molecule
     */
    void regenerateCoordinates();

    /**
     * Undoably clear the molecule and non-molecular objects.
     */
    void clear();

    /**
     * Undoably select or deselect the specified atoms, bonds, and non-molecular
     * objects.
     *
     * @param atoms The atoms to select or deselect
     * @param bonds The bonds to select or deselect
     * @param non_molecular_objects The non-molecular objects to select or
     * deselect
     * @param select_mode Whether to select, deselect, toggle selection, or
     * select-only (i.e. clear the selection and then select)
     *
     * @note If either an attachment point dummy atom or the associated bond are
     * selected (or deselected, etc.), then the other will automatically be
     * selected (or deselected, etc.) as well.
     */
    void select(const std::unordered_set<const RDKit::Atom*>& atoms,
                const std::unordered_set<const RDKit::Bond*>& bonds,
                const std::unordered_set<const NonMolecularObject*>&
                    non_molecular_objects,
                const SelectMode select_mode);

    /**
     * Undoably clear all selected atoms, bonds, and non-molecular objects.
     */
    void clearSelection();

    /**
     * Select all atoms, bonds, and non-molecular objects.
     */
    void selectAll();

    /**
     * Select all unselected objects, and deselect all selected objects.
     */
    void invertSelection();

    /**
     * Remove all selected atoms, bonds, and non-molecular objects
     */
    void removeSelected();

    /**
     * Mutate a set of atoms atoms to the specified element.
     * @param atoms The atoms to mutate
     */
    void mutateAtoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                     const Element element);

    /**
     * Mutate a set of atoms to the specified query.
     */
    void mutateAtoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                     const AtomQuery atom_query);

    /**
     * Mutate a set of atoms to the specified atom.
     */
    void mutateAtoms(const std::unordered_set<const RDKit::Atom*>& from_atoms,
                     const RDKit::Atom& to_atom);

    /**
     * Change the R-group of existing atoms. This method can also mutate
     * non-R-group atoms into R-group atoms. If an explicit number is not
     * specified, the next available R-group number will be used.
     */
    void mutateRGroups(const std::unordered_set<const RDKit::Atom*>& atoms);
    void mutateRGroups(const std::unordered_set<const RDKit::Atom*>& atoms,
                       const unsigned int r_group_num);

    /**
     * Mutate all selected bonds
     * @param bond_tool The type of bond to mutate to
     */
    void mutateBonds(const std::unordered_set<const RDKit::Bond*>& bonds,
                     BondTool bond_tool);

    /**
     * Adjust the charge on all target atoms by the specified value. To
     * decrement the charge, pass in a negative number.
     */
    void
    adjustChargeOnAtoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                        const int increment_by);

    /**
     * Adjust the number of unpaired electrons on all target atoms by the
     * specified value. To decrement the charge, pass in a negative number.
     */
    void adjustRadicalElectronsOnAtoms(
        const std::unordered_set<const RDKit::Atom*>& atoms, int increment_by);

    /**
     * If any atoms of the set have implicit hydrogens, add explicit hydrogens
     * to all atoms of the set. Otherwise, hide all explicit hydrogens on all
     * atoms of the set.
     * @param atoms The atoms to consider.
     */
    void toggleExplicitHsOnAtoms(
        const std::unordered_set<const RDKit::Atom*>& atoms);

  signals:

    /**
     * Signal emitted when the molecule is changed (e.g. by adding or removing
     * atoms or bonds), excluding selection changes, which will trigger
     * selectionChanged instead and coordinates changes which will trigger
     * coordinatesChanged.
     */
    void moleculeChanged();

    /**
     * Signal emitted when the non-molecular objects are changed, excluding
     * selection changes, which will trigger selectionChanged instead, and
     * coordinate changes, which will trigger coordinatesChanged.
     */
    void nonMolecularObjectsChanged();

    /**
     * Signal emitted when the coordinates of the molecule or the non-molecular
     * objects are changed
     */
    void coordinatesChanged();

    /**
     * Signal emitted when selection is changed.  Note that atoms, bonds, and
     * non-molecular objects will automatically be deselected when they are
     * removed and reselected if the removal is undone, and this signal *will*
     * be emitted in those scenario.
     */
    void selectionChanged();

    /**
     * Signal emitted when a reaction arrow is added.  Note that this signal
     * will be emitted *in addition to* nonMolecularObjectsChanged.
     */
    void reactionArrowAdded();

  protected:
    RDKit::RWMol m_mol = RDKit::RWMol();
    std::vector<NonMolecularObject> m_pluses;
    std::optional<NonMolecularObject> m_arrow;
    std::unordered_map<int, NonMolecularObject*> m_tag_to_non_molecular_object;

    int m_next_atom_tag = 0;
    int m_next_bond_tag = 0;
    int m_next_non_molecular_tag = 0;
    std::unordered_set<int> m_selected_atom_tags;
    std::unordered_set<int> m_selected_bond_tags;
    std::unordered_set<int> m_selected_non_molecular_tags;

    std::vector<HighlightingInfo> m_highlighting_info;

    /**
     * create an empty conformer for m_mol so it's ready to be used by other
     * functions. This is called in the constructor and whenever the model is
     * cleared.
     */
    void initializeMol();

    /**
     * Transform all coordinates using the given function.
     * @param desc The description to use for the redo/undo command.
     * @param function The function to use to transform the coordinates.
     * @param merge_id The merge id to use for the redo/undo command. If this is
     * different from -1, the command will be merged with the previous command
     * if they share the same merge id.
     * @param atoms The atoms to transform. If both atoms and
     * non_molecular_objects are empty, then all objects will be transformed.
     * @param non_molecular_objects The non-molecular objects to transform. If
     * both atoms and non_molecular_objects are empty, then all objects will be
     * transformed.
     */
    void transformCoordinatesWithFunction(
        const QString& desc, std::function<void(RDGeom::Point3D&)> function,
        MergeId merge_id = MergeId::NO_MERGE,
        std::unordered_set<const RDKit::Atom*> atoms = {},
        std::unordered_set<const NonMolecularObject*> non_molecular_objects =
            {});

    void addExplicitHs(const std::unordered_set<const RDKit::Atom*>& atoms);
    void removeExplicitHs(const std::unordered_set<const RDKit::Atom*>& atoms);

    /**
     * Set the atom tag for the specified atom
     */
    void setTagForAtom(RDKit::Atom* const atom, const int atom_tag);

    /**
     * Find the atom tag for the specified atom.  The passed in value must not
     * be nullptr.
     * @param atom The atom to get the tag for
     * @param allow_null If this is true, return -1 if the nullptr is passed in.
     * Otherwise, raise an exception
     *
     * Note that this method is const, but not marked as such because RDKit
     * does not implement a const version of getAtomBookmarks.  (And because
     * this method isn't public, so it's not worth using const_cast unless we
     * need this to be const.)
     */
    int getTagForAtom(const RDKit::Atom* const atom,
                      const bool allow_null = false);

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
     * Helpers that return a non-const pointer to the specified atom or bond
     * by round-tripping through bookmarks of the model's underlying mol.
     *
     * Because MolModel's public APIs return const pointers to m_mol and its
     * contents, the scene and all widgets only have access to these non-const
     * pointers. When they are then passed back to MolModel (regardless of
     * whether the method wants to change the underlying mol or not), we need
     * to accept const pointers, confirm that it is in fact part of the
     * underlying mol, and work on non-const pointers.
     */
    RDKit::Atom* getMutableAtom(const RDKit::Atom* const atom);
    RDKit::Bond* getMutableBond(const RDKit::Bond* const bond);

    /**
     * Get the atomic number and element name for the specified element.
     */
    std::pair<unsigned int, QString> getAddElementInfo(const Element& element);

    /**
     * Undoably select or deselect the specified atoms and bonds.
     *
     * @param atom_tags Tags for the atoms to select or deselect
     * @param bond_tags Tags for the bonds to select or deselect
     * @param non_molecular_tags Tags for the non-molecular objects to select or
     * deselect
     * @param select_mode Whether to select, deselect, toggle selection, or
     * select-only (i.e. clear the selection and then select)
     *
     * @note This method does *not* expand attachment points (i.e. make sure
     * that both the dummy atom and the bond are included).  That is instead
     * handled in the select method.
     */
    void selectTags(const std::unordered_set<int>& atom_tags,
                    const std::unordered_set<int>& bond_tags,
                    const std::unordered_set<int>& non_molecular_tags,
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
     * @return The atom, bond, and non-molecular tags for all objects that are
     * not currently selected
     */
    std::tuple<std::unordered_set<int>, std::unordered_set<int>,
               std::unordered_set<int>>
    getAllUnselectedTags();

    /**
     * Undoably select or deselect the specified atoms and bonds.  If
     * selecting, all given tags must be currently unselected.  If deselecting,
     * all given tags must be currently selected.
     *
     * @param filtered_atom_tags Tags for the atoms to select or deselect.
     * @param filtered_bond_tags Tags for the bonds to select or deselect.
     * @param filtered_non_molecular_tags Tags for the non-molecular objects to
     * select or deselect.
     * @param to_select Whether to select or deselect the specified atoms and
     * bonds.
     * @param description The description for the undo command
     */
    void doSelectionCommand(
        const std::unordered_set<int>& filtered_atom_tags,
        const std::unordered_set<int>& filtered_bond_tags,
        const std::unordered_set<int>& filtered_non_molecular_tags,
        const bool to_select, const QString& description);

    /**
     * Execute the given function immediately, and use snapshots to provide undo
     * and redo functionality for the changes that the function carries out.
     * Note that doing, undoing, or redoing this command will invalidate all
     * RDKit Atom and Bond objects.
     *
     * @param do_func The function to call for the initial do.  This function
     * should not emit any signals, as that will be handled by the snapshots.
     * @param description A description for the undo command
     * @param to_be_changed Whether do_func updates the RDKit molecule,
     * the non-molecular objects, or both.  If this value includes
     * WhatChanged::MOLECULE, then update_molecule_on_change() will be called
     * automatically after do_func is executed.  If this value includes
     * WhatChanged::NON_MOL_OBJS, then updateNonMolecularMetadata() will be
     * called automatically.
     */
    void doCommandUsingSnapshots(const std::function<void()> do_func,
                                 const QString& description,
                                 const WhatChangedType to_be_changed);

    /**
     * @return a copy of the MolModel state
     */
    MolModelSnapshot takeSnapshot() const;

    /**
     * Update the current MolModel state to match the provided snapshot
     * @param snapshot The snapshot to restore
     * @param what_changed Whether this snapshot updates the RDKit molecule,
     * the non-molecular objects, or both
     * @param selection_changed Whether this snapshot changes the selection
     * @param arrow_added Whether this snapshot change adds a reaction arrow.
     */
    void restoreSnapshot(const MolModelSnapshot& snapshot,
                         const WhatChangedType what_changed,
                         const bool selection_changed, const bool arrow_added);

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
     * atoms in the chain.  This method must only be called as part of an undo
     * command.
     *
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
    void addAtomChainCommandFunc(const AtomFunc create_atom,
                                 const std::vector<RDGeom::Point3D>& coords,
                                 const BondFunc create_bond,
                                 const int bound_to_atom_tag);

    /**
     * Add a non-molecular object (a plus sign or a reaction arrow).  This
     * method must only be called as part of an undo command.
     * @param type The type of non-molecular object to add
     * @param coords The coordinates for the new non-molecular object
     */
    void addNonMolecularObjectCommandFunc(const NonMolecularType& type,
                                          const RDGeom::Point3D& coords);

    /**
     * Update m_tag_to_non_molecular_object so that it reflects the currently
     * existant non-molecular objects
     */
    void updateNonMolecularMetadata();

    /**
     * Remove an atom from the molecule.  This method must only be called as
     * part of an undo command.
     * @param atom_tag The tag of the atom to delete
     * @return whether selection was changed by this action
     */
    bool removeAtomCommandFunc(const int atom_tag);

    /**
     * Add a bond to the molecule.  This method must only be called as part of
     * an undo command.
     * @param start_atom_tag The tag of the bond's start atom
     * @param end_atom_tag The tag of the bond's end atom
     * @param create_bond A function that will return the new bond to add.  Note
     * that this method will add a *copy* of the returned bond, as the
     * shared_ptr is responsible for deleting the returned bond.
     */
    void addBondCommandFunc(const int start_atom_tag, const int end_atom_tag,
                            const BondFunc create_bond);

    /**
     * Remove a bond from the molecule.  This method must only be called as part
     * of an undo command.
     * @return whether selection was changed by this action
     */
    bool removeBondCommandFunc(const int bond_tag, const int start_atom_tag,
                               const int end_atom_tag);

    /**
     * Remove a non-molecular object from the model.  This method must only be
     * called as part of an undo command.
     * @param cur_tag The tag of the non-molecular object to delete
     * @return whether selection was changed by this action
     */
    bool removeNonMolecularObjectCommandFunc(const int cur_tag);

    /**
     * set new coordinates for a set of atoms and non-molecular objects.  This
     * method must only be called as part of an undo command.
     */
    void setCoordinates(
        const std::vector<const RDKit::Atom*>& atoms,
        const std::vector<RDGeom::Point3D>& atom_coords,
        const std::vector<const NonMolecularObject*> non_molecular_objects,
        const std::vector<RDGeom::Point3D>& non_mol_coords);

    /**
     * Remove the specified atoms and bonds from the molecule.  This method
     * must only be called as part of an undo command.
     * @param atom_tags The atom tags to delete
     * @param bond_tags_with_atoms A list of (bond tag, bond's start atom tag,
     * bond's end atom tag) for the bonds to delete
     * @param sgroups The substance groups to delete
     * @param non_molecular_tags The non-molecular tags to delete
     */
    void removeCommandFunc(
        const std::vector<int>& atom_tags,
        const std::vector<std::tuple<int, int, int>>& bond_tags_with_atoms,
        const std::unordered_set<const RDKit::SubstanceGroup*>& sgroups,
        const std::vector<int>& non_molecular_tags);

    /**
     * Add all atoms and bonds from the given molecule into this molecule. This
     * method must only be called as part of an undo command.
     *
     * @param mol The molecule to add
     */
    void addMolCommandFunc(const RDKit::ROMol& mol);

    /**
     * Add the given reaction to the model. This method must only be called as
     * part of an undo command.
     */
    void addReactionCommandFunc(const RDKit::ChemicalReaction& reaction);

    /**
     * Change the element of an existing atom, or change a query atom to a
     * non-query atom.  This method must only be called as part of an undo
     * command.
     * @param atom_tag The tag of the atom to mutate
     * @param create_atom A function that will return a new atom to mutate to.
     * Note that this method will add a *copy* of the returned atom, as the
     * shared_ptr is responsible for deleting the returned atom.
     */
    void mutateAtomCommandFunc(const int atom_tag, const AtomFunc create_atom);

    /**
     * Change the type of an existing bond, or change a query bond to a
     * non-query bond.  This method must only be called as part of an undo
     * command.
     * @param bond_tag The tag of the bond to mutate
     * @param create_bond A function that will return a new bond to mutate to.
     * Note that this method will add a *copy* of the returned bond, as the
     * shared_ptr is responsible for deleting the returned bond.
     */
    void mutateBondCommandFunc(const int bond_tag, const BondFunc create_bond);

    /**
     * Set the mapping number of a set of atoms.  This method must only be
     * called from an undo command.
     * @param atom_tags The tag of the atoms to add mappings to
     * @param mapping The new mapping number of the atom. If this is 0 the
     * mapping is removed
     */
    void setAtomMappingCommandFunc(const std::unordered_set<int>& atom_tags,
                                   const int mapping);

    /**
     * Add explicit hydrogens to an existing atom. This method must only be
     * called as part of an undo command.
     */
    void addExplicitHsCommandFunc(
        const std::unordered_set<const RDKit::Atom*>& atom_tags);

    /**
     * Remove explicit hydrogens from an existing atom.  This method must only
     * be called from an undo command.
     */
    void removeExplicitHsCommandFunc(
        const std::unordered_set<const RDKit::Atom*>& atom_tags);

    /**
     * Reverse an existing bond (i.e. swap the start and end atoms).  This
     * method must only be called as part of an undo command.
     */
    void flipBondCommandFunc(const RDKit::Bond* bond);

    /**
     * Clear the molecule.  This method must only be called as part of an undo
     * command.
     */
    void clearCommandFunc();

    /**
     * Select or deselect the specified atoms and bonds.  This method must only
     * be called as part of an undo command.
     *
     * @param atom_tags The atom tags to select or deselect
     * @param bond_tags The bond tags to select or deselect
     * @param non_molecular_tags The non-molecular tags to select or deselect
     * @param selected Whether to select or deselect the specified atoms and
     * bonds
     */
    void
    setSelectionCommandFunc(const std::unordered_set<int>& atom_tags,
                            const std::unordered_set<int>& bond_tags,
                            const std::unordered_set<int>& non_molecular_tags,
                            const bool selected);

    /**
     * Clear all selected atoms and bonds.  This method must only be called
     * as part of an undo command.
     */
    void clearSelectionCommandFunc();

    /**
     * Add a fragment to the molecule and create new bonds between the new
     * fragment and the core (i.e. the existing molecule).  This method must
     * only be called as part of an undo command.
     * @param fragment The fragment to add
     * @param mutations_to_core_atoms A list of core atoms to mutate, where each
     *     atom is given as (core atom index, function that returns a bond to
     *     mutate to)
     * @param mutations_to_core_bonds A list of core bonds to mutate, where each
     * bond is given as (core bond index, function that returns a bond to mutate
     * to)
     * @param additions_to_core_bonds A list of bonds to make between two core
     * atoms, where each bond is given as a tuple of
     *   - atom index for the starting atom
     *   - atom index for the ending atom
     *   - function that returns a bond instance for the new bond
     * @param core_to_frag_bonds_by_idx Bonds to make between the core and the
     * fragment, formatted as a map of {core atom index: {fragment atom index:
     * info about bond to create}}, and bond info is given as a tuple of
     *   - function that returns a bond instance for the new bond
     *   - whether the bond starts with the fragment atom (as opposed to
     *     starting with the core atom)
     */
    void addFragmentCommandFunc(
        const RDKit::ROMol& fragment,
        const std::vector<std::pair<unsigned int, AtomFunc>>&
            mutations_to_core_atoms,
        const std::vector<std::pair<unsigned int, BondFunc>>&
            mutations_to_core_bonds,
        const std::vector<std::tuple<unsigned int, unsigned int, BondFunc>>&
            additions_to_core_bonds,
        const AtomIdxToFragBondMap& core_to_frag_bonds_by_idx);

    /**
     * Mutate a set of atoms to an Atom object returned by the provided
     * function.
     * @param atoms The atoms to mutate
     */
    void mutateAtoms(const std::unordered_set<const RDKit::Atom*>& atoms,
                     const AtomFunc& create_atom);

    /**
     * Mutate a set of bonds to a Bond object returned by the provided
     * function.
     * @param bonds The bonds to mutate
     * @param flip_matching_bonds If true, then any bonds that already match
     * both bond_type and bond_dir will be flipped (i.e. their start and end
     * atoms will be switched).  If false, then any bonds that already match
     * will be ignored.
     */
    void mutateBonds(const std::unordered_set<const RDKit::Bond*>& bonds,
                     const BondFunc& create_bond,
                     bool flip_matching_bonds = false);

    /**
     * Combine pairs of atoms into single atoms with all bonds from both
     * original atoms.  This is typically used to merge overlapping atoms at the
     * end of a rotation or translation.  This method must only be called as
     * part of an undo command.
     *
     * @param atoms_to_merge A vector of {atom to keep, stationary atom} pairs.
     * The resulting merged atom has the coordinates of the stationary atom, but
     * all other properties (e.g. element, charge) are taken from the atom to
     * keep.
     */
    void mergeAtomsCommandFunc(
        const std::vector<std::pair<const RDKit::Atom*, const RDKit::Atom*>>&
            atoms_to_merge);

    /**
     * Clean up coordinates for the current reaction.  This method must only be
     * called when an arrow is present in the workspace, and as part of an undo
     * command.
     */
    void regenerateReactionCoordinatesCommandFunc();

    /**
     * Create an RDKit reaction object containing this model's contents.  This
     * method must only be called when an arrow is present in the workspace.
     *
     * @param strip_tags Whether to remove all internal MolModel atom and bond
     * properties (e.g. atom and bond tags).  These properties should always be
     * removed before passing the reaction to the user, but may be helpful if
     * the reaction object is for internal MolModel use only.
     */
    boost::shared_ptr<RDKit::ChemicalReaction>
    createReaction(const bool strip_tags) const;
};

} // namespace sketcher
} // namespace schrodinger
