#pragma once

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>

namespace schrodinger
{
namespace rdkit_extensions
{

const std::string CHAIN_HYDROGEN{"chainHydrogen"};
const std::string SGROUP_PROP_CLASS{"CLASS"};
const std::string SGROUP_PROP_LABEL{"LABEL"};
const std::string SGROUP_PROP_RESNUM{"ID"};
const std::string SGROUP_PROP_TYPE{"TYPE"};

/**
 * Translate SUP (Super Atom) substance groups into PDB residue
 * info on atoms so that atomistic -> monomeric conversion
 * can be done.
 *
 * Since there is no guarentee that SUP groups cover the entire structure,
 * atoms that are not part of a SUP group will placed into separate monomers
 * by connectivity. Hydrogen leaving groups will be removed and not translated
 * into a separate monomer. Connectivity information is NOT determined by the
 * attachment point labels on the SUP groups, but by the actual atomistic
 * connectivity.
 *
 * @param mol The molecule to process
 * @return true if the molecule had SUP groups and was processed
 */
bool processSupGroups(RDKit::ROMol& mol);

/**
 * Create SUP substance groups from PDB residue info on a molecule's atoms
 *
 * This function creates SUP (superatom) substance groups based on the
 * PDB residue information stored on each atom's MonomerInfo. Atoms with
 * the same chain ID and residue number are grouped together.
 *
 * These SUP groups can be used to preserve residue information when writing
 * to SDF.
 *
 * @param mol molecule to add SUP groups to (modified in place)
 */
void createSupGroupsFromResidueInfo(RDKit::RWMol& mol);

} // namespace rdkit_extensions
} // namespace schrodinger
