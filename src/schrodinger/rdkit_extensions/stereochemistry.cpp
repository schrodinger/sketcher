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

std::string get_atom_chirality_label(const RDKit::Atom& atom, bool strip_abs)
{
    // SKETCH-1729: non-CIP ranked atoms have no priority; do not show any label
    auto has_non_CIP_neighbor = [](const auto& atom) {
        auto stereo_info = RDKit::Chirality::detail::getStereoInfo(&atom);
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
    if (has_non_CIP_neighbor(atom)) {
        return "";
    }

    auto possible_but_not_specified = [](const auto& atom) {
        int possible = 0;
        atom.getPropIfPresent(RDKit::common_properties::_ChiralityPossible,
                              possible);
        return !atom.hasProp(RDKit::common_properties::_CIPCode) && possible;
    };
    if (possible_but_not_specified(atom)) {
        return "(?)";
    }

    std::string chirality;
    atom.getPropIfPresent<std::string>(RDKit::common_properties::atomNote,
                                       chirality);

    if (strip_abs && chirality.find(ABSOLUTE_STEREO_PREFIX) == 0) {
        // remove the ABSOLUTE_STEREO_PREFIX
        chirality = chirality.substr(ABSOLUTE_STEREO_PREFIX.size());
    }
    return chirality;
}

std::string get_simplified_stereo_annotation(const RDKit::ROMol& mol)
{
    auto sgs = mol.getStereoGroups();
    if (sgs.size() != 1)
        return "";
    boost::dynamic_bitset<> chiralAts(mol.getNumAtoms());
    for (const auto atom : mol.atoms()) {
        if (atom->getChiralTag() > RDKit::Atom::ChiralType::CHI_UNSPECIFIED &&
            atom->getChiralTag() < RDKit::Atom::ChiralType::CHI_OTHER) {
            chiralAts.set(atom->getIdx(), 1);
        }
    }
    for (const auto atm : sgs[0].getAtoms()) {
        chiralAts.set(atm->getIdx(), 0);
    }
    if (!chiralAts.any())
        return "";
    // all specified chiral centers are accounted for by this StereoGroup.
    if (sgs[0].getGroupType() == RDKit::StereoGroupType::STEREO_OR) {
        return "OR enantiomer";
    } else if (sgs[0].getGroupType() == RDKit::StereoGroupType::STEREO_AND) {
        return "AND enantiomer";
    }
    return "";
}

std::string get_bond_stereo_label(const RDKit::Bond& bond)
{
    std::string label;
    bond.getPropIfPresent(RDKit::common_properties::bondNote, label);
    return label;
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