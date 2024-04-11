#pragma once

#include "schrodinger/rdkit_extensions/definitions.h"

#include <boost/shared_ptr.hpp>
#include <string>
#include <string_view>
#include <vector>

const std::string ANNOTATION{"ANNOTATION"};
const std::string LINKAGE{"attachmentPoints"};
const std::string ATOM_LABEL{"atomLabel"};
const std::string SUPPLEMENTARY_INFORMATION{"SUPPLEMENTARY_INFORMATION"};
const std::string BRANCH_LINKAGE{"R3-R1"};
const std::string BACKBONE_LINKAGE{"R2-R1"};
const std::string HELM_MODEL{"HELM_MODEL"};
const std::string MONOMER_LIST{"MONOMER_LIST"};
const std::string UNKNOWN_MONOMER{"UNKNOWN_MONOMER"};

// NOTE: These are to allow replacement of the python api
const std::string BRANCH_MONOMER{"isBranchMonomer"};
const std::string SMILES_MONOMER{"isSmilesMonomer"};
const std::string CUSTOM_BOND{"customBond"};
const std::string EXTENDED_ANNOTATIONS{"extended_annotations"};

const std::string ARM_PAIR_KEY{"Antibody arms"};
const std::string STRAND_PAIR_KEY{"Double Helix"};

namespace RDKit
{
class ROMol;
class RWMol;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{
[[nodiscard]] RDKIT_EXTENSIONS_API bool
is_coarse_grain_mol(const RDKit::ROMol& mol);

// Helper apis to get atom indices belonging to the inputs polymers
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<unsigned int>
get_atoms_in_polymer_chain(const RDKit::ROMol& mol,
                           std::string_view polymer_id);

// NOTE: Duplicate polymer ids will be ignored
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<unsigned int>
get_atoms_in_polymer_chains(const RDKit::ROMol& mol,
                            const std::vector<std::string_view>& polymer_ids);

// Extracts polymers with the provided ids from a coarse grain mol. Polymer
// information that span across non-selected polymers will be deleted.
//
// NOTE: Polymer groups are unsupported
[[nodiscard]] RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
extract_helm_polymers(const RDKit::ROMol& mol,
                      const std::vector<std::string_view>& polymer_ids);

} // namespace rdkit_extensions
} // namespace schrodinger
