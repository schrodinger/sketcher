#include "schrodinger/rdkit_extensions/stereochemistry.h"

#include <rdkit/CIPLabeler/CIPLabeler.h>
#include <rdkit/CIPLabeler/TooManyNodesException.h>
#include <rdkit/GraphMol/Chirality.h>
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
    RDKit::MolOps::assignStereochemistry(mol);
}

void assign_CIP_labels(RDKit::RWMol& mol)
{
    try {
        // This number of calculation cycles takes:
        // ~1s on a Linux Intel(R) Xeon(R) W-2123 CPU @ 3.60GHz.
        // ~0.5s on a 2019 Mac Book Pro with a Intel i7 @ 2.6 GHz.
        // ~1.5s using the WASM sketcher on either of these.
        unsigned max_cycles = 1000000;

        RDKit::CIPLabeler::assignCIPLabels(mol, max_cycles);
    } catch (const RDKit::CIPLabeler::MaxIterationsExceeded&) {
        // CIP label calculation "timed out". Some labels will be omitted.
    } catch (const RDKit::CIPLabeler::TooManyNodesException&) {
        // CIP label calculation graph became too big. It's unlikely we hit
        // this, we'll most probably hit the "time out' before. Still,
        // keep any labels we found so far.
    }
}

std::string get_atom_chirality_label(const RDKit::Atom& atom)
{
    std::string chirality;
    auto stereo_info = RDKit::Chirality::detail::getStereoInfo(&atom);

    // non-CIP ranked atoms have no priority, so potentially chiral
    // atoms should always show an 'undefined' label (SKETCH-1825)
    auto has_non_CIP_neighbor = [](const auto& atom, const auto& stereo_info) {
        auto& mol = atom.getOwningMol();
        for (auto atom_idx : stereo_info.controllingAtoms) {
            auto controlling_atom = mol.getAtomWithIdx(atom_idx);
            if (RDKit::getAtomRLabel(controlling_atom) != 0 ||
                schrodinger::rdkit_extensions::is_attachment_point_dummy(
                    *controlling_atom)) {
                return true;
            }
        }
        return false;
    };

    if (has_non_CIP_neighbor(atom, stereo_info)) {
        return ""; // SKETCH-1729: Don't show a defined label
    }

    if (int possible = 0;
        !atom.getPropIfPresent(RDKit::common_properties::_CIPCode, chirality) &&
        atom.getPropIfPresent(RDKit::common_properties::_ChiralityPossible,
                              possible) &&
        possible) {
        return "?"; // possible, but not specified
    }
    return chirality;
}

std::string get_bond_stereo_label(const RDKit::Bond& bond)
{
    std::string label;
    // bond stereo can never be undefined because we always
    // have 2d coords, so double bonds will always be either
    // STEREOANY (no label), E/Z or non-stereo capable (no label).
    bond.getPropIfPresent(RDKit::common_properties::_CIPCode, label);
    return label;
}

void add_enhanced_stereo_to_chiral_atoms(RDKit::ROMol& mol)
{
    // translateChiralFlagToStereoGroups() works based on the chiral flag, so
    // make sure it is set. We default to chiral flag on (absolute stereo) if
    // the flag is not present, or else we'll flip stereo on things that don't
    // come from .sdf.
    int chiral_flag{1};
    if (!mol.getPropIfPresent(RDKit::common_properties::_MolFileChiralFlag,
                              chiral_flag)) {
        mol.setProp(RDKit::common_properties::_MolFileChiralFlag, chiral_flag);
    }

    RDKit::translateChiralFlagToStereoGroups(mol);
}

void reapply_molblock_wedging(RDKit::ROMol& rdk_mol)
{
    for (auto& bond : rdk_mol.bonds()) {
        // only change the wedging if the bond was wedged in the input data (we
        // recognize that by looking for the properties the mol file parser
        // sets):
        unsigned int bond_dir_val{0};
        if (bond->getPropIfPresent(RDKit::common_properties::_MolFileBondStereo,
                                   bond_dir_val)) {
            // v2000
            if (bond_dir_val == 1) {
                bond->setBondDir(RDKit::Bond::BondDir::BEGINWEDGE);
            } else if (bond_dir_val == 6) {
                bond->setBondDir(RDKit::Bond::BondDir::BEGINDASH);
            } else if (bond_dir_val == 4) {
                bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
            }
        } else if (bond->getPropIfPresent(
                       RDKit::common_properties::_MolFileBondCfg,
                       bond_dir_val)) {
            // v3000
            if (bond_dir_val == 1) {
                bond->setBondDir(RDKit::Bond::BondDir::BEGINWEDGE);
            } else if (bond_dir_val == 3) {
                bond->setBondDir(RDKit::Bond::BondDir::BEGINDASH);
            } else if (bond_dir_val == 2) {
                bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
            }
        } else {
            auto bond_dir = bond->getBondDir();
            if (bond_dir == RDKit::Bond::BondDir::BEGINWEDGE ||
                bond_dir == RDKit::Bond::BondDir::BEGINDASH ||
                bond_dir == RDKit::Bond::BondDir::UNKNOWN) {
                bond->setBondDir(RDKit::Bond::BondDir::UNKNOWN);
            }
        }
    }
}

void wedgeMolBonds(RDKit::ROMol& mol, const RDKit::Conformer* conf)
{
    std::vector<RDKit::Bond*> attachment_dummy_bonds;
    for (auto atom : mol.atoms()) {
        if (rdkit_extensions::is_attachment_point_dummy(*atom)) {
            auto bonds = mol.atomBonds(atom);
            auto bond = *bonds.begin();
            if (bond->getBondType() == RDKit::Bond::SINGLE) {
                attachment_dummy_bonds.push_back(bond);
                bond->setBondType(RDKit::Bond::OTHER);
            }
        }
    }

    RDKit::ClearSingleBondDirFlags(mol);
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