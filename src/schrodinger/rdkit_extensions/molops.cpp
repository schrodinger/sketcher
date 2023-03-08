/* -------------------------------------------------------------------------
 * Implements schrodinger::rdkit_extensions:: miscellaneous mol operations
 mol conversion
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */
#include "schrodinger/rdkit_extensions/molops.h"

#include <GraphMol/MolOps.h>

#include "schrodinger/rdkit_extensions/constants.h"

namespace schrodinger
{
namespace rdkit_extensions
{

void apply_sanitization(RDKit::RWMol& mol, Sanitization sanitization)
{
    RDLog::LogStateSetter silence_rdkit_logging;

    unsigned failed_op{0};
    int ops = RDKit::MolOps::SANITIZE_ALL;
    if (sanitization == Sanitization::PARTIAL) {
        ops = RDKit::MolOps::SANITIZE_SYMMRINGS |
              RDKit::MolOps::SANITIZE_SETCONJUGATION |
              RDKit::MolOps::SANITIZE_SETHYBRIDIZATION |
              RDKit::MolOps::SANITIZE_ADJUSTHS;
    }
    RDKit::MolOps::sanitizeMol(mol, failed_op, ops);
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

bool is_attachment_point_dummy(const RDKit::Atom& atom)
{
    std::string label;
    return atom.getAtomicNum() == 0 && atom.getTotalDegree() == 1 &&
           atom.getPropIfPresent(RDKit::common_properties::atomLabel, label) &&
           label.find(ATTACHMENT_POINT_LABEL_PREFIX) == 0;
}

void add_enhanced_stereo_to_chiral_atoms(RDKit::ROMol& mol)
{
    auto stereo_groups = mol.getStereoGroups();
    std::unordered_set<RDKit::Atom*> seen_chiral_atoms;
    for (const auto& sg : stereo_groups) {
        const auto& atoms = sg.getAtoms();
        seen_chiral_atoms.insert(atoms.begin(), atoms.end());
    }

    std::vector<RDKit::Atom*> ungrouped_atoms;
    for (auto atom : mol.atoms()) {
        if ((atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CW ||
             atom->getChiralTag() == RDKit::Atom::CHI_TETRAHEDRAL_CCW) &&
            seen_chiral_atoms.count(atom) == 0) {
            ungrouped_atoms.push_back(atom);
        }
    }

    if (ungrouped_atoms.empty()) {
        return;
    }

    // Default to chiral flag on (absolute stereo) if the flag is not present,
    // or else we'll flip stereo on things that don't come from .sdf.
    // For .sdf, we explicitly set chiral_flag=0 if the flag is not present
    // (see convert.cpp). Note that only adds the flag property when it is on,
    // when it is off, no property is present.
    int chiral_flag{1};
    mol.getPropIfPresent(RDKit::common_properties::_MolFileChiralFlag,
                         chiral_flag);

    if (chiral_flag == 1) {
        // We have a chiral flag: append to or create the ABS group
        // (ABS is unique!)

        auto is_abs = [](const auto& sg) {
            return sg.getGroupType() == RDKit::StereoGroupType::STEREO_ABSOLUTE;
        };
        auto abs_group =
            std::find_if(stereo_groups.begin(), stereo_groups.end(), is_abs);

        if (abs_group != stereo_groups.end()) {
            auto abs_atoms = abs_group->getAtoms();
            ungrouped_atoms.insert(ungrouped_atoms.end(), abs_atoms.begin(),
                                   abs_atoms.end());
            stereo_groups.erase(abs_group);
        }

        stereo_groups.emplace_back(RDKit::StereoGroupType::STEREO_ABSOLUTE,
                                   std::move(ungrouped_atoms));
    } else {
        // No chiral flag: ungrouped atoms go into a new AND group
        stereo_groups.emplace_back(RDKit::StereoGroupType::STEREO_AND,
                                   std::move(ungrouped_atoms));
    }

    mol.setStereoGroups(stereo_groups);
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

void removeHs(RDKit::RWMol& rdk_mol)
{
    RDKit::MolOps::RemoveHsParameters ps;

    // We always remove H on queries; for sketcher import, all atoms are created
    // as QueryAtoms as that they might be changed into queries later on; for
    // conversion from 3D Structure, there is no way to create queries.
    ps.removeWithQuery = true;

    // Disable displaying warnings
    ps.showWarnings = false;

    bool sanitize = false;
    RDKit::MolOps::removeHs(rdk_mol, ps, sanitize);
}

} // namespace rdkit_extensions
} // namespace schrodinger
