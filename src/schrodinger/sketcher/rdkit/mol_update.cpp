#include "schrodinger/sketcher/rdkit/mol_update.h"

#include <rdkit/GraphMol/Atropisomers.h>
#include <rdkit/GraphMol/CIPLabeler/CIPLabeler.h>
#include <rdkit/GraphMol/CIPLabeler/TooManyNodesException.h>
#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/FileParsers/MolFileStereochem.h>
#include <rdkit/GraphMol/MolOps.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/sketcher/rdkit/coord_utils.h"
#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/sketcher/rdkit/sgroup.h"
#include "schrodinger/rdkit_extensions/stereochemistry.h"
#include "schrodinger/rdkit_extensions/constants.h"
#include "schrodinger/sketcher/rdkit/atom_properties.h"

#include <boost/format.hpp>

namespace schrodinger
{
namespace sketcher
{

/**
 * Sketcher-specific assignment of enhanced stereo to all chiral centers. This
 * preserves any assigned enhanced stereo group and their numbering, while
 * assigning new enhanced stereo groups to any chiral centers that lack them.
 */
static void add_enhanced_stereo_to_chiral_atoms(RDKit::ROMol& mol)
{
    // If the chiral flag property is present (i.e. MDL), leave it as-is;
    // otherwise, force the chiral flag "on" (i.e. SMILES, MAE, etc.)
    // The chiral flag state controls how enhanced stereo is assigned below,
    // where the "off" state assigns AND and the "on" state assigns ABS to
    // ungrouped chiral centers.
    int chiral_flag{1};
    if (!mol.getPropIfPresent(RDKit::common_properties::_MolFileChiralFlag,
                              chiral_flag)) {
        mol.setProp(RDKit::common_properties::_MolFileChiralFlag, chiral_flag);
    }

    // Assign enhanced stereo groups to all centers that lack it
    // NOTE: this clears the chiral flag property on the mol
    RDKit::translateChiralFlagToStereoGroups(mol);
}

/**
 * @return whether the given mol has CFG/stereo properties present that would
 * have been read from MDL molblock input
 */
static bool has_molblock_cfgs(RDKit::ROMol& mol)
{
    auto bonds = mol.bonds();
    return std::any_of(bonds.begin(), bonds.end(), [](auto& bond) {
        return bond->getBondType() == RDKit::Bond::SINGLE &&
               (bond->hasProp(RDKit::common_properties::_MolFileBondCfg) ||
                bond->hasProp(RDKit::common_properties::_MolFileBondStereo));
    });
}

/**
 * Workaround for RDKit's assignChiralTypesFromBondDirs setting explicit H
 * counts on 3-coordinate chiral atoms (see Chirality.cpp:3533).
 * This wrapper preserves the original explicit H counts and restores them
 * after the call.
 * TODO: After updating with rdkit/rdkit#8872, investigate removing this
 */
static void
assign_chiral_types_from_bond_dirs_preserving_explicit_hs(RDKit::RWMol& mol)
{
    std::vector<unsigned int> initial_num_explicit_hs;
    initial_num_explicit_hs.reserve(mol.getNumAtoms());
    for (auto atom : mol.atoms()) {
        initial_num_explicit_hs.push_back(atom->getNumExplicitHs());
    }

    RDKit::MolOps::assignChiralTypesFromBondDirs(mol);

    // Reset explicit H counts back to their original values and recalculate
    // implicit H counts.
    for (size_t i = 0; auto atom : mol.atoms()) {
        unsigned int initial_explicit_hs = initial_num_explicit_hs[i++];
        if (atom->getNumExplicitHs() != initial_explicit_hs) {
            atom->setNumExplicitHs(initial_explicit_hs);
            atom->updatePropertyCache(/* strict = */ false);
        }
    }
}

/**
 * Uses coordinates and bond directions to calculate stereochemistry from
 * scratch.
 */
static void
assign_stereochemistry_with_bond_directions_and_coordinates(RDKit::RWMol& mol)
{
    // Reset the calculated intermediate stereo markers RDKit uses in its stereo
    // assignment so we can recalculate them from scratch. Also store
    // UP/DOWN/WIGGLY bond directions so we can restore them after cleanIt=true.
    for (auto atom : mol.atoms()) {
        atom->setChiralTag(RDKit::Atom::CHI_UNSPECIFIED);
    }
    std::vector<RDKit::Bond::BondDir> bond_dirs;
    bond_dirs.reserve(mol.getNumBonds());
    for (auto bond : mol.bonds()) {
        bond_dirs.push_back(bond->getBondDir());
        if (bond->getBondDir() == RDKit::Bond::BondDir::ENDDOWNRIGHT ||
            bond->getBondDir() == RDKit::Bond::BondDir::ENDUPRIGHT) {
            bond->setBondDir(RDKit::Bond::BondDir::NONE);
        }

        // assignStereochemistry(cleanIt=true) will not clear
        // atropisomer stereo markers, so we need to do it manually.
        bond->setStereo(RDKit::Bond::BondStereo::STEREONONE);
    }

    // Convert up/down bonds into parities
    assign_chiral_types_from_bond_dirs_preserving_explicit_hs(mol);

    // Detect atropisomer chirality & clean up unwanted/redundant atrop markers:
    // This is currently not done by assignStereochemistry()
    auto conf = mol.getConformer();
    RDKit::Atropisomers::detectAtropisomerChirality(mol, &conf);
    RDKit::MolOps::cleanupAtropisomers(mol);

    // Calculate stereo bond cues from coordinates. RDKit can't calculate
    // stereo bonds without these.
    RDKit::MolOps::setDoubleBondNeighborDirections(mol, &conf);

    // Assign stereo from those parities (primarily for bond stereo)
    rdkit_extensions::UseModernStereoPerception use_modern_stereo;

    bool cleanIt = true;
    bool force = true;
    bool flagPossibleStereoCenters = true;
    RDKit::MolOps::assignStereochemistry(mol, cleanIt, force,
                                         flagPossibleStereoCenters);

    // Restore bond directions. We want to preserve EITHERDOUBLE since
    // it's what we use to draw crossed bonds.
    auto bond_it = bond_dirs.begin();
    for (auto bond : mol.bonds()) {
        if (auto bond_dir = *bond_it++;
            bond_dir == RDKit::Bond::BondDir::BEGINDASH ||
            bond_dir == RDKit::Bond::BondDir::BEGINWEDGE ||
            bond_dir == RDKit::Bond::BondDir::EITHERDOUBLE ||
            bond_dir == RDKit::Bond::BondDir::UNKNOWN) {
            bond->setBondDir(bond_dir);
        }
    }
}

/**
 * Sketcher-specific assignment of CIP labels, which bypasses the assignment
 * after a certain number of cycles or certain exceptions are hit. For this
 * function to work properly, the input mol must have had stereo assigned.
 */
static void assign_CIP_labels(RDKit::RWMol& mol)
{
    try {
        // Set in SHARED-11140
        unsigned max_cycles = 2000000;

        RDKit::CIPLabeler::assignCIPLabels(mol, max_cycles);
    } catch (const RDKit::CIPLabeler::MaxIterationsExceeded&) {
        // CIP label calculation "timed out". Some labels will be omitted.
    } catch (const RDKit::CIPLabeler::TooManyNodesException&) {
        // CIP label calculation graph became too big. It's unlikely we hit
        // this, we'll most probably hit the "time out' before. Still,
        // keep any labels we found so far.
    }
}

// If we can't resolve an atoms CIP code (e.g. due to timing out),
// the {cip} placeholder in the label won't be updated to the actual
// code, and the placeholder will be displayed. To prevent this from
// happening, we need to postprocess the labels and remove the ones that
// haven't been replaced.
// We only care about ABS groups because these are the only ones for
// which we show the CIP code in the label (AND/OR groups just show the
// type of the group, which we always know, and for ungrouped atoms we
// only add a label if the CIP code is known).
static void clear_abs_labels_with_unresolved_cip(RDKit::RWMol& mol)
{
    for (const auto& sg : mol.getStereoGroups()) {
        // We only care about the ABS group
        if (sg.getGroupType() == RDKit::StereoGroupType::STEREO_ABSOLUTE) {
            for (auto atom : sg.getAtoms()) {
                if (!atom->hasProp(RDKit::common_properties::_CIPCode)) {
                    atom->clearProp(RDKit::common_properties::atomNote);
                }
            }
        }
    }
}

/**
 * Called once when a mol is first brought into the sketcher
 */
void prepare_mol(RDKit::ROMol& mol)
{
    // Make sure valences are available.
    mol.updatePropertyCache(false);

    // Add 2D coordinates only if the molecule does not already have them
    // present (ie specified via molblock, SMILES extension, etc.) Note that
    // this call will intentionally keep a 3d conformer if one is found, since
    // Maestro may use that later to align the Sketcher structure with the
    // Maestro workspace structure
    update_2d_coordinates(mol);

    // Update all input chiral centers to have enhanced stereo; honors MDL input
    add_enhanced_stereo_to_chiral_atoms(mol);

    // Convert parities back to wedges/dashes
    if (has_molblock_cfgs(mol)) {
        RDKit::Chirality::reapplyMolBlockWedging(mol);
    } else { // Otherwise recalculate chiral bond directions
        rdkit_extensions::wedgeMolBonds(mol, &mol.getConformer());
    }
}

/**
 * Called every time the mol is changed in the sketcher
 */
void update_molecule_on_change(RDKit::RWMol& mol)
{
    // save the simplified stereo annotation as a mol property
    RDKit::Chirality::simplifyEnhancedStereo(
        mol, /*removeAffectedStereoGroups = */ false);

    // Explicitly update the brackets for the sgroups
    update_s_group_brackets(mol);

    // DO NOT REPLACE THIS WITH rdkit_extensions::apply_sanitization(PARTIAL) or
    // any sanitization that includes FINDRADICALS: it will turn atoms with
    // unsatisfied valences (coming from SMARTS) into RADICALs, and we don't
    // want that. See SKETCH-2190 or test_smarts_no_radicals in test_mol_model
    // for more info.
    bool no_strict = false;
    mol.updatePropertyCache(no_strict);
    // update hybridization states so that coordinate calculations produce the
    // correct geometry
    RDKit::MolOps::setConjugation(mol);
    RDKit::MolOps::setHybridization(mol);
    // This is the bare minimum we need to do to be able to detect stereo.
    RDKit::MolOps::symmetrizeSSSR(mol);

    assign_stereochemistry_with_bond_directions_and_coordinates(mol);

    // Generate R/S/E/Z labels
    assign_CIP_labels(mol);

    // Add the newly found chiral atoms to StereoGroups
    add_enhanced_stereo_to_chiral_atoms(mol);

    // Preserve IDs of any new enhanced stereo groups, which includes any
    // groups newly inserted from file/paste input. Run this after
    // cleanupStereoGroups(), because that one only preserves input ids
    RDKit::forwardStereoGroupIds(mol);

    /**
     * note that the ABSOLUTE_STEREO_PREFIX might be stripped away by the
     * rendering code depending on rendering preferences. Note that the molModel
     * remains unaware of these rendering preferences (which are part of
     * sketcherModel) because we want everything in the molModel to be savable
     * as an undoable snapshot. The rendering preferences might get exposed to
     * the GUI in the future and we don't want them to be undoable.
     */
    std::string abs_label =
        rdkit_extensions::ABSOLUTE_STEREO_PREFIX + "({cip})";
    std::string or_label = rdkit_extensions::OR_STEREO_PREFIX + "{id}";
    std::string and_label = rdkit_extensions::AND_STEREO_PREFIX + "{id}";

    // addStereoAnnotations will not clear existing annotations,
    // just override the ones on atoms/bonds that require labels.
    // So we need to clear them ourselves first to get dir of
    // outdated labels.
    for (auto atom : mol.atoms()) {
        atom->clearProp(RDKit::common_properties::atomNote);
    }
    for (auto bond : mol.bonds()) {
        bond->clearProp(RDKit::common_properties::bondNote);
    }

    RDKit::Chirality::addStereoAnnotations(mol, abs_label, or_label, and_label);

    clear_abs_labels_with_unresolved_cip(mol);
}

} // namespace sketcher
} // namespace schrodinger
