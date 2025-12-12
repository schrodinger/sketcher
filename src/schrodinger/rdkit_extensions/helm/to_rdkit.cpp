#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/range/irange.hpp>
#include <charconv>
#include <fmt/format.h>
#include <iterator>
#include <memory>
#include <string>
#include <vector>

#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/Atom.h>
#include <rdkit/GraphMol/QueryAtom.h>
#include <rdkit/GraphMol/Bond.h>
#include <rdkit/GraphMol/SubstanceGroup.h>

#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/helm_parser.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"

using namespace schrodinger::rdkit_extensions;

/*
 * The structure of the output romol is as follows:
 *      * polymers
 *          * are represented as COP sgroups with included atoms and bonds
 *          * annotations are stored as the ANNOTATION property on the sgroup
 *      * monomers:
 *          * are represented as atoms
 *          * monomer ids are stored as the ATOM_LABEL property on each atom
 *          * annotations are stored as the ANNOTATION property on an annotated
 *            atom
 *          * repeated monomers with ambiguous repetitions are stored as SRU
 *          sgroups
 *          * NOTE: A hash of the monomer id is stored as the isotope number
 *            for the monomer. This allow us to use RDKit substructure
 *            search apis for monomeric romols.
 *      * linkages:
 *          * are represented as bonds between atoms
 *          * annotations are stored as the ANNOTATION property on an annotated
 *            connection
 *          * the participating rgroups are stored as the LINKAGE property on
 *            each bond
 *          * hydrogen bonds are depicted as zero-order bonds and other types
 *            of linkages as single bonds
 *
 *      * polymer groups and extended annotations
 *          * are store as a DAT sgroup with the SUPPLEMENTARY_INFORMATION field
 *            name
 *          * the sgroup is a vector with two strings encoding the polymer
 *            groups and extended annotations
 *          * these are currently lumped together since there's no obvious use
 *            for them other than roundtripping to a HELMV2.0 string.
 *
 * TODO: Add example use cases
 */

namespace
{

template <class PropType, class RDKitObject> PropType get_property
    [[nodiscard]] (const RDKitObject* obj, const std::string& propname)
{
    return obj->template getProp<PropType>(propname);
}
} // namespace

namespace helm
{

namespace
{
[[nodiscard]] std::vector<::RDKit::RWMol>
create_helm_polymers(const std::vector<polymer>& polymers);

[[nodiscard]] ::RDKit::RWMol create_helm_polymer(const polymer& polymer);

void add_helm_connections(::RDKit::RWMol& mol,
                          const std::vector<connection>& connections);

void add_polymer_groups_and_extended_annotations(
    ::RDKit::RWMol& mol, const std::vector<polymer_group> polymer_groups,
    const std::string_view& extended_annotations);

} // namespace

[[nodiscard]] std::unique_ptr<::RDKit::RWMol>
helm_to_rdkit(const std::string& helm_string, bool do_throw)
{
    auto parsed_info = helm::parse_helm(helm_string, do_throw);
    if (!parsed_info) {
        return nullptr;
    }

    const auto polymers = create_helm_polymers(parsed_info->polymers);

    std::unique_ptr<::RDKit::RWMol> mol(new ::RDKit::RWMol);
    std::for_each(polymers.begin(), polymers.end(),
                  [&](const auto& polymer) { mol->insertMol(polymer); });

    // TODO: Keep track of unused rgroups
    add_helm_connections(*mol, parsed_info->connections);

    add_polymer_groups_and_extended_annotations(
        *mol, parsed_info->polymer_groups, parsed_info->extended_annotations);

    // store this as an easy way to identify monomeric mols
    mol->setProp<bool>(HELM_MODEL, true);
    return mol;
}

namespace
{
[[nodiscard]] std::vector<::RDKit::RWMol>
create_helm_polymers(const std::vector<helm::polymer>& polymers)
{
    std::vector<::RDKit::RWMol> polymer_mols;
    polymer_mols.reserve(polymers.size());
    std::transform(polymers.begin(), polymers.end(),
                   std::back_inserter(polymer_mols), create_helm_polymer);
    return polymer_mols;
}

void condense_monomer_list(const std::string_view& polymer_id,
                           std::string_view& monomer_id)
{
    // monomer lists with wildcard or unknown monomers can be condensed into
    // a single monomer
    static constexpr std::string_view wildcard_monomer{"*"};
    static constexpr std::string_view unknown_amino_acid{"X"};
    static constexpr std::string_view unknown_nucleotide{"N"};
    static constexpr unsigned char peptide_prefix{'P'};
    static constexpr unsigned char nucleotide_prefix{'R'};

    const auto wildcard_pos = monomer_id.find(wildcard_monomer);
    if (wildcard_pos != std::string_view::npos) {
        monomer_id = wildcard_monomer;
        return;
    }
    // We prevent missing monomers in the parser
    std::string_view unknown_monomer_id{"_"};
    if (polymer_id.front() == peptide_prefix) {
        unknown_monomer_id = unknown_amino_acid;
    } else if (polymer_id.front() == nucleotide_prefix) {
        unknown_monomer_id = unknown_nucleotide;
    }

    const auto unknown_monomer_pos = monomer_id.find(unknown_monomer_id);
    if (unknown_monomer_pos != std::string_view::npos) {
        monomer_id = unknown_monomer_id;
        return;
    }
    // NOTE: if we can't condense it, then just display the first residue name
    monomer_id = monomer_id.substr(0, monomer_id.find_first_of(",+"));
    if (monomer_id.front() == '[') {
        monomer_id = monomer_id.substr(1, monomer_id.size() - 2);
    }
}

[[nodiscard]] std::string get_monomer_id(const std::string_view& polymer_id,
                                         const monomer& monomer)
{
    // Strip square brackets from multi-character monomers and try to condense
    // monomer lists. I'm not doing this in separate clauses because the entry
    // [dA],G,[dC] could result in unintentionally stripping the square brackets
    auto monomer_id = monomer.id;
    if (monomer.is_list) {
        condense_monomer_list(polymer_id, monomer_id);
    }
    // Inline SMILES tokens are stripped in the parser.
    else if (!monomer.is_smiles && monomer_id.front() == '[' &&
             monomer_id.back() == ']') {
        monomer_id = monomer.id.substr(1, monomer.id.size() - 2);
    }
    return std::string{monomer_id};
}

[[nodiscard]] ::RDKit::RWMol create_helm_polymer(const polymer& polymer)
{
    // NOTE: Not making this static because gcc complains it's unused
    auto create_ambiguous_repetition_sru = [](auto& mol, const auto& repetition,
                                              std::string_view polymer_id) {
        ::RDKit::SubstanceGroup sru_sgroup(&mol, "SRU");

        const auto& [start, size, num_repetitions, annotation] = repetition;
        const auto repeated_atoms = boost::irange(start, start + size);
        sru_sgroup.setAtoms({repeated_atoms.begin(), repeated_atoms.end()});

        // NOTE: SRU sgroups only allow 2 or 4 bonds so we should handle this
        unsigned int bond1 = 0; // prev bond to non-repeated group
        unsigned int bond2 = 0; // next bond to non-repeated group
        // Do the end first, so we don't have to care
        // about the start adding a QueryAtom or not
        if (start + size == mol.getNumAtoms()) {
            auto idx = mol.addAtom(new ::RDKit::QueryAtom(), true, true);
            mol.getAtomWithIdx(idx)->setProp(REPETITION_DUMMY_ID,
                                             std::string{polymer_id});
            addConnection(mol, start + size - 1, idx, BACKBONE_LINKAGE);
            bond2 = mol.getNumBonds() - 1;
        } else {
            bond2 = start + size - 1;
        }
        if (start == 0) {
            auto idx = mol.addAtom(new ::RDKit::QueryAtom(), true, true);
            mol.getAtomWithIdx(idx)->setProp(REPETITION_DUMMY_ID,
                                             std::string{polymer_id});
            addConnection(mol, 0, idx, BACKBONE_LINKAGE);
            bond1 = mol.getNumBonds() - 1;
        } else {
            bond1 = start - 1;
        }
        sru_sgroup.setBonds({bond1, bond2});

        // NOTE: Doing this since enumeration requires an HH or HT
        // CONNECTION
        sru_sgroup.setProp("CONNECT", "HT");
        // NOTE: Set number of repetitions as LABEL to allow recognition by
        // enumeration apis
        sru_sgroup.setProp("LABEL", std::string{num_repetitions});
        if (!annotation.empty()) {
            sru_sgroup.setProp(ANNOTATION, std::string{annotation});
        }
        return sru_sgroup;
    };

    std::vector<::RDKit::SubstanceGroup> monomer_lists;
    ::RDKit::RWMol mol;
    int prev_backbone_idx = -1;
    unsigned int residue_number = 1;
    for (const auto& monomer : polymer.monomers) {
        auto monomer_type =
            monomer.is_smiles ? MonomerType::SMILES : MonomerType::REGULAR;
        auto idx = addMonomer(mol, get_monomer_id(polymer.id, monomer),
                              residue_number, polymer.id, monomer_type);
        ++residue_number;
        auto atom = mol.getAtomWithIdx(idx);
        if (!monomer.annotation.empty()) {
            atom->setProp(ANNOTATION, std::string{monomer.annotation});
        }
        // NOTE: Doing it this way since it has no functional use at the
        // moment.
        if (monomer.is_list) {
            atom->setProp(MONOMER_LIST, std::string{monomer.id});
        }

        if (prev_backbone_idx != -1) {
            const auto connection_type = monomer.is_branch
                                             ? ConnectionType::SIDECHAIN
                                             : ConnectionType::FORWARD;
            addConnection(mol, prev_backbone_idx, idx, connection_type);
        }
        prev_backbone_idx = (monomer.is_branch ? prev_backbone_idx : idx);
    }

    // if the monomer sequence has ambiguous repetitions, create an SRU
    // sgroup to store that information
    std::vector<::RDKit::SubstanceGroup> repeated_monomers;
    const auto& repetitions = polymer.repetitions;
    std::transform(
        repetitions.begin(), repetitions.end(),
        std::back_inserter(repeated_monomers), [&](const auto& repetition) {
            return create_ambiguous_repetition_sru(mol, repetition, polymer.id);
        });

    if (!polymer.annotation.empty()) {
        RDKit::SubstanceGroup annotation_sgroup{&mol, "COP"};
        annotation_sgroup.setProp("ID", std::string{polymer.id});
        annotation_sgroup.setProp(ANNOTATION, std::string{polymer.annotation});
        ::RDKit::addSubstanceGroup(mol, annotation_sgroup);
    }

    std::for_each(repeated_monomers.begin(), repeated_monomers.end(),
                  [&](const auto& sru_sgroup) {
                      ::RDKit::addSubstanceGroup(mol, sru_sgroup);
                  });
    // monomer lists
    std::for_each(
        monomer_lists.begin(), monomer_lists.end(),
        [&](const auto& sgroup) { ::RDKit::addSubstanceGroup(mol, sgroup); });

    return mol;
}

// <from_atom_idx, to_atom_idx, linkage_property, annotation>
using linkage_info =
    std::tuple<unsigned int, unsigned int, std::string, std::string>;

std::vector<linkage_info>
enumerate_connection_monomers(const ::RDKit::RWMol& mol,
                              const connection& connection)
{

    // Utility function to map residue numbers to atom indices
    static auto get_monomers =
        [](const auto& mol, const auto& polymer_id,
           const auto& residues) -> std::vector<unsigned int> {
        const auto polymer = get_polymer(mol, polymer_id);

        std::vector<std::string> residue_list;
        boost::split(residue_list,
                     std::string{residues.begin(), residues.end()},
                     boost::is_any_of(",+"));

        std::vector<unsigned int> atom_indices;
        for (const auto& residue : residue_list) {
            auto is_residue_number =
                std::all_of(residue.begin(), residue.end(),
                            [&](unsigned char c) { return std::isdigit(c); });
            if (is_residue_number) {
                atom_indices.push_back(polymer.atoms[std::stoi(residue) - 1]);
            } else {
                // if it's a residue name, we have to look for all atoms with
                // said residue name. Select everything for wildcard values
                std::copy_if(polymer.atoms.begin(), polymer.atoms.end(),
                             std::back_inserter(atom_indices),
                             [&](const auto& idx) {
                                 auto* atom = mol.getAtomWithIdx(idx);
                                 return (residue.find_first_of("*?") !=
                                             std::string::npos ||
                                         get_property<std::string>(
                                             atom, ATOM_LABEL) == residue);
                             });
            }
        }
        return atom_indices;
    };

    const auto source_monomers =
        get_monomers(mol, connection.from_id, connection.from_res);
    const auto target_monomers =
        get_monomers(mol, connection.to_id, connection.to_res);

    std::vector<linkage_info> connected_monomers;
    const std::string linkage =
        fmt::format("{}-{}", connection.from_rgroup, connection.to_rgroup);
    std::for_each(source_monomers.begin(), source_monomers.end(),
                  [&](const auto source_monomer) {
                      std::for_each(
                          target_monomers.begin(), target_monomers.end(),
                          [&](const auto target_monomer) {
                              connected_monomers.push_back(
                                  {source_monomer, target_monomer, linkage,
                                   std::string{connection.annotation}});
                          });
                  });
    return connected_monomers;
}

void add_helm_connections(::RDKit::RWMol& mol,
                          const std::vector<connection>& connections)
{
    std::vector<linkage_info> connected_monomers;
    std::for_each(
        connections.begin(), connections.end(), [&](const auto& connection) {
            auto bonds_info = enumerate_connection_monomers(mol, connection);
            connected_monomers.insert(connected_monomers.end(),
                                      bonds_info.begin(), bonds_info.end());
        });

    std::ranges::for_each(connected_monomers, [&](const auto& info) {
        const auto& [from_atom, to_atom, linkage, annotation] = info;
        constexpr bool is_custom_bond = true;
        addConnection(mol, from_atom, to_atom, linkage, is_custom_bond);
        if (!annotation.empty()) {
            mol.getBondBetweenAtoms(from_atom, to_atom)
                ->setProp(ANNOTATION, annotation);
        }
    });
}

void add_polymer_groups_and_extended_annotations(
    ::RDKit::RWMol& mol, const std::vector<polymer_group> polymer_groups,
    const std::string_view& extended_annotations)
{
    // nothing to do here
    if (polymer_groups.empty() && extended_annotations.empty()) {
        return;
    }

    ::RDKit::SubstanceGroup sgroup(&mol, "DAT");
    sgroup.setProp("FIELDNAME", SUPPLEMENTARY_INFORMATION);

    // convert polymer groups to HELMV2.0-style strings
    std::vector<std::string> polymer_group_strings;
    std::transform(polymer_groups.begin(), polymer_groups.end(),
                   std::back_inserter(polymer_group_strings),
                   [](const auto& polymer_group) {
                       return fmt::format("{}({})", polymer_group.id,
                                          polymer_group.items);
                   });

    std::vector<std::string> datafields{
        boost::algorithm::join(polymer_group_strings, "|"),
        std::string{extended_annotations}};
    sgroup.setProp("DATAFIELDS", datafields);
    ::RDKit::addSubstanceGroup(mol, sgroup);

    // NOTE: For python apis
    if (!extended_annotations.empty()) {
        mol.setProp(EXTENDED_ANNOTATIONS, std::string(extended_annotations));
    }
}
} // namespace
} // namespace helm
