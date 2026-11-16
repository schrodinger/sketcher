#include "schrodinger/rdkit_extensions/stereochemistry.h"

#include <rdkit/GraphMol/GraphMol.h>
#include <rdkit/GraphMol/Chirality.h>
#include <rdkit/GraphMol/StereoGroup.h>
#include <rdkit/GraphMol/FileParsers/MolFileStereochem.h>

#include "schrodinger/rdkit_extensions/molops.h"
#include "schrodinger/rdkit_extensions/rgroup.h"

namespace schrodinger
{
namespace rdkit_extensions
{

UseModernStereoPerception::UseModernStereoPerception() :
    m_stereo_algo_state{RDKit::Chirality::getUseLegacyStereoPerception()}
{
    RDKit::Chirality::setUseLegacyStereoPerception(false);
}

UseModernStereoPerception::~UseModernStereoPerception()
{
    RDKit::Chirality::setUseLegacyStereoPerception(m_stereo_algo_state);
}

void assign_stereochemistry(RDKit::ROMol& mol)
{
    if (mol.getNumConformers() != 0) {
        auto conf = mol.getConformer();
        const auto& ps = conf.getPositions();

        if (std::any_of(ps.begin(), ps.end(),
                        [](auto p) { return p.length() > 1e-4; })) {
            // If we have all-zero coordinates, we don't care about
            // the 2D/3D flag (which is read, rather than checked)

            if (conf.is3D()) {
                // If we have 3D coordinates, calculate both chirality
                // and stereo bonds from them.
                RDKit::MolOps::assignStereochemistryFrom3D(mol);
                return;

            } else {
                // If we have 2D coordinates available, use them
                // to refresh stereo bond directions
                RDKit::MolOps::setDoubleBondNeighborDirections(mol, &conf);
            }
        }
    }

    // The general case, valid both no conformation and calculate
    // stereochemistry from parity, and stereo bonds from bond directions.
    bool cleanIt = false;
    bool force = true;
    bool flagPossibleStereoCenters = true;
    RDKit::MolOps::assignStereochemistry(mol, cleanIt, force,
                                         flagPossibleStereoCenters);
}

bool is_wedge_or_dash(const RDKit::Bond& bond)
{
    auto bond_dir = bond.getBondDir();
    return bond_dir == RDKit::Bond::BondDir::BEGINWEDGE ||
           bond_dir == RDKit::Bond::BondDir::BEGINDASH;
}

/**
 * @return whether the given atom has a tetrahedral chiral tag, which is the
 * only type of chirality that RDKit's WedgeMolBonds() adds wedges for
 */
static bool has_tetrahedral_chiral_tag(const RDKit::Atom& atom)
{
    auto chiral_tag = atom.getChiralTag();
    return chiral_tag == RDKit::Atom::CHI_TETRAHEDRAL_CW ||
           chiral_tag == RDKit::Atom::CHI_TETRAHEDRAL_CCW;
}

/**
 * Update the existing bond directions of the given molecule so that they can be
 * kept when the molecule is passed to RDKit's WedgeMolBonds(). See the
 * wedgeMolBonds() docstring for a description of the bond directions that are
 * changed.
 */
static void clean_up_existing_bond_directions(RDKit::ROMol& mol,
                                              const RDKit::Conformer* conf)
{
    std::vector<bool> atom_has_wedge(mol.getNumAtoms(), false);
    for (auto bond : mol.bonds()) {
        if (is_wedge_or_dash(*bond)) {
            atom_has_wedge[bond->getBeginAtomIdx()] = true;
        }
    }

    for (auto bond : mol.bonds()) {
        auto bond_dir = bond->getBondDir();
        if (bond->getBondType() == RDKit::Bond::SINGLE &&
            (bond_dir == RDKit::Bond::BondDir::ENDUPRIGHT ||
             bond_dir == RDKit::Bond::BondDir::ENDDOWNRIGHT)) {
            // these only describe the stereochemistry of a neighboring double
            // bond, which is recorded on the double bond itself
            bond->setBondDir(RDKit::Bond::BondDir::NONE);
        } else if (!is_wedge_or_dash(*bond)) {
            continue;
        } else if (has_tetrahedral_chiral_tag(*bond->getBeginAtom())) {
            // recalculate the direction of this bond from the chiral tag, which
            // will flip the bond if it disagreed with the tag. (Note that
            // wedgeBond leaves the bond as-is if it's not a single bond or if
            // it can't determine a direction.)
            RDKit::Chirality::wedgeBond(bond, bond->getBeginAtomIdx(), conf);
        } else if (has_tetrahedral_chiral_tag(*bond->getEndAtom()) &&
                   !atom_has_wedge[bond->getEndAtomIdx()]) {
            // RDKit won't add a wedge to a chiral center if a wedge or dash
            // ends there, even though such a bond says nothing about the
            // center's stereochemistry
            bond->setBondDir(RDKit::Bond::BondDir::NONE);
        }
    }
}

/**
 * Remove the MDL wedging info for wedged and dashed single bonds. The info for
 * wavy single bonds and crossed double bonds is left in place, since RDKit's
 * CXSMILES and MAE writers read it.
 */
static void clear_molblock_wedge_and_dash_info(RDKit::ROMol& mol)
{
    auto clear_prop_if_wedge_or_dash = [](RDKit::Bond* bond, const auto& prop,
                                          const int wedge_value,
                                          const int dash_value) {
        int value = 0;
        if (bond->getPropIfPresent(prop, value) &&
            (value == wedge_value || value == dash_value)) {
            bond->clearProp(prop);
        }
    };
    for (auto bond : mol.bonds()) {
        if (bond->getBondType() != RDKit::Bond::SINGLE) {
            continue;
        }
        clear_prop_if_wedge_or_dash(
            bond, RDKit::common_properties::_MolFileBondCfg, 1, 3);
        clear_prop_if_wedge_or_dash(
            bond, RDKit::common_properties::_MolFileBondStereo, 1, 6);
    }
}

void wedgeMolBonds(RDKit::ROMol& mol, const RDKit::Conformer* conf,
                   bool keep_existing_wedges)
{
    std::vector<RDKit::Bond*> attachment_dummy_bonds;
    for (auto atom : mol.atoms()) {
        if (is_attachment_point_dummy(*atom)) {
            auto bonds = mol.atomBonds(atom);
            auto bond = *bonds.begin();
            if (bond->getBondType() == RDKit::Bond::SINGLE) {
                attachment_dummy_bonds.push_back(bond);
                bond->setBondType(RDKit::Bond::OTHER);
            }
        }
    }

    if (keep_existing_wedges) {
        clean_up_existing_bond_directions(mol, conf);
        // The wedges and dashes are stored as bond directions, so discard any
        // MDL wedging info for them. Otherwise, that info would go stale as
        // soon as a wedge is edited, and it would be copied along with the
        // molecule (e.g. the Sketcher's clipboard pickle includes it).
        clear_molblock_wedge_and_dash_info(mol);
    } else {
        RDKit::ClearSingleBondDirFlags(mol);
        RDKit::Chirality::clearMolBlockWedgingInfo(mol);
    }
    // Note that RDKit's WedgeMolBonds() won't add a wedge to any chiral atom
    // that already has one

    try {
        // Temporarily silence RDKit's loggers
        RDLog::LogStateSetter silence_rdkit_logging;

        RDKit::WedgeMolBonds(mol, conf);
    } catch (const Invar::Invariant&) {
        // RDkit wasn't able to find a 'wedgeable' bond for some chiral atom
    }

    // Restore the dummies
    for (auto bond : attachment_dummy_bonds) {
        bond->setBondType(RDKit::Bond::SINGLE);
    }
}

} // namespace rdkit_extensions
} // namespace schrodinger
