#include "schrodinger/rdkit_extensions/stereochemistry.h"

#include <rdkit/GraphMol/GraphMol.h>
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
    bool cleanIt = false;
    bool force = true;
    bool flagPossibleStereoCenters = true;
    RDKit::MolOps::assignStereochemistry(mol, cleanIt, force,
                                         flagPossibleStereoCenters);
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
                is_attachment_point_dummy(*controlling_atom)) {
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

    // TODO: make sure new groups have IDS
}

void wedgeMolBonds(RDKit::ROMol& mol, const RDKit::Conformer* conf)
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