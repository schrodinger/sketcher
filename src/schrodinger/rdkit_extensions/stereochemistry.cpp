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
    RDKit::Chirality::clearMolBlockWedgingInfo(mol);

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
