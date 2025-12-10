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
const std::string REPETITION_DUMMY_ID{"REPETITION_DUMMY_ID"};

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
class Atom;
class Bond;
class SubstanceGroup;
} // namespace RDKit

namespace schrodinger
{
namespace rdkit_extensions
{

// acts as a view to a polymer chain, and is only valid as long as the owning
// mol is valid.
struct Chain {
    std::vector<unsigned int> atoms;
    std::vector<unsigned int> bonds;
    std::string annotation;
    // std::string polymer_id;
};

[[nodiscard]] RDKIT_EXTENSIONS_API bool isMonomeric(const RDKit::ROMol& mol);

// Helper apis to get atom indices belonging to the inputs polymers
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<unsigned int>
get_atoms_in_polymer_chain(const RDKit::ROMol& mol,
                           std::string_view polymer_id);

// NOTE: Duplicate polymer ids will be ignored
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<unsigned int>
get_atoms_in_polymer_chains(const RDKit::ROMol& mol,
                            const std::vector<std::string_view>& polymer_ids);

// Extracts polymers with the provided ids from a monomer mol. Polymer
// information that span across non-selected polymers will be deleted.
//
// NOTE: Polymer groups are unsupported
[[nodiscard]] RDKIT_EXTENSIONS_API boost::shared_ptr<RDKit::RWMol>
extract_helm_polymers(const RDKit::ROMol& mol,
                      const std::vector<std::string_view>& polymer_ids);

[[nodiscard]] RDKIT_EXTENSIONS_API std::string
get_polymer_id(const ::RDKit::Atom* atom);

[[nodiscard]] RDKIT_EXTENSIONS_API unsigned int
get_residue_number(const ::RDKit::Atom* atom);

[[nodiscard]] RDKIT_EXTENSIONS_API bool
is_dummy_atom(const ::RDKit::Atom* atom);

/**
 * Connections (in the HELM context) are bonds that are written in the
 * connection section of a HELM string. These are bonds that either connect
 * two polymers together or complete a loop in a polymer.
 *
 * In order to call this function, chain and residue information needs to be
 * assigned to the atoms. If this CG molecule was not created using the HELM
 * parser, assign_chains needs to be called first.
 *
 * @param monomer_mol monomeric molecule to extract connections from
 * @return All 'connection' bonds that will be in the connection section of
 * a HELM string
 */
[[nodiscard]] RDKIT_EXTENSIONS_API std::vector<unsigned int>
get_connections(const ::RDKit::ROMol& monomer_mol);

std::vector<std::string>
    RDKIT_EXTENSIONS_API get_polymer_ids(const RDKit::ROMol& monomer_mol);

Chain RDKIT_EXTENSIONS_API get_polymer(const RDKit::ROMol& monomer_mol,
                                       std::string_view polymer_id);

/**
 * Get the HELM supplementary info substance group
 *
 * @param cg_mol coarse grain molecule
 * @return a pointer to the supplementary info substance group if one exists
 */
[[nodiscard]] RDKIT_EXTENSIONS_API const RDKit::SubstanceGroup*
get_supplementary_info(const RDKit::ROMol& cg_mol);

/**
 * @return whether the given substance group contains HELM supplementary info
 */
[[nodiscard]] RDKIT_EXTENSIONS_API bool
is_supplementary_information_s_group(const RDKit::SubstanceGroup& sgroup);

/**
 * @return whether the given substance group contains a polymer annotation
 */
[[nodiscard]] RDKIT_EXTENSIONS_API bool
is_polymer_annotation_s_group(const RDKit::SubstanceGroup& sgroup);

} // namespace rdkit_extensions
} // namespace schrodinger
