#include "schrodinger/rdkit_extensions/fasta/to_string.h"

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <fmt/format.h>
#include <fmt/ranges.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <string_view>

#include "schrodinger/rdkit_extensions/fasta/monomers.h"
#include "schrodinger/rdkit_extensions/helm.h"

using schrodinger::rdkit_extensions::get_connections;
using schrodinger::rdkit_extensions::get_polymer;
using schrodinger::rdkit_extensions::get_polymer_id;
using schrodinger::rdkit_extensions::get_polymer_ids;
using schrodinger::rdkit_extensions::get_residue_number;

static const std::string SGROUP_TYPE{"TYPE"};

static void check_for_unsupported_features(const ::RDKit::ROMol& mol);

static void
check_nucleotide_structure(const RDKit::ROMol& mol,
                           const schrodinger::rdkit_extensions::Chain& chain);

[[nodiscard]] static std::string
get_peptide_fasta(const RDKit::ROMol& mol,
                  const schrodinger::rdkit_extensions::Chain& chain);

[[nodiscard]] static std::string
get_nucleotide_fasta(const RDKit::ROMol& mol,
                     const schrodinger::rdkit_extensions::Chain& chain);

template <class PropType, class RDKitObject> [[nodiscard]] static auto
get_property(const RDKitObject* obj, const std::string& propname)
{
    return obj->template getProp<PropType>(propname);
}

[[nodiscard]] static std::string get_monomer_id(const ::RDKit::Atom* atom)
{
    if (atom->hasProp(MONOMER_LIST)) {
        // we do this sorting since we hard code supported monomer information.
        // this allow us to identify (A+G) and (G+A) as the same monomer. For
        // the purposes of FASTA conversion, we are not distinguishing between
        // unions and mutually exclusive lists i.e. (A+G) and (A,G)
        std::vector<std::string> monomer_list;
        boost::split(monomer_list,
                     get_property<std::string>(atom, MONOMER_LIST),
                     boost::is_any_of("+,"));
        std::sort(monomer_list.begin(), monomer_list.end());
        return fmt::format("({})", fmt::join(monomer_list, "+"));
    } else if (atom->getProp<bool>(SMILES_MONOMER)) {
        return "X";

    } else {
        const auto id = get_property<std::string>(atom, ATOM_LABEL);
        return (id.size() == 1 ? id : fmt::format("[{}]", id));
    }
}

[[nodiscard]] static bool is_biopolymer(std::string_view polymer_id)
{
    // PEPTIDE[0-9]* or RNA[0-9]*
    return polymer_id.front() == 'P' || polymer_id.front() == 'R';
}

[[nodiscard]] static std::string
get_fasta_from_biopolymer(const RDKit::ROMol& mol,
                          const schrodinger::rdkit_extensions::Chain& chain,
                          std::string_view polymer_id)
{
    // PEPTIDE[0-9]*
    if (polymer_id.front() == 'P') {
        return get_peptide_fasta(mol, chain);
    }
    return get_nucleotide_fasta(mol, chain);
}

namespace fasta
{

[[nodiscard]] std::string rdkit_to_fasta(const ::RDKit::ROMol& mol)
{
    check_for_unsupported_features(mol);

    fmt::memory_buffer output_fasta;
    // constexpy bool sanitize = false;
    for (const auto& polymer_id : get_polymer_ids(mol)) {
        auto chain = get_polymer(mol, polymer_id);

        fmt::format_to(std::back_inserter(output_fasta), "{}\n",
                       get_fasta_from_biopolymer(mol, chain, polymer_id));
    }

    return {output_fasta.data(), output_fasta.size()};
}

} // namespace fasta

static void check_for_unsupported_features(const ::RDKit::ROMol& mol)
{
    if (!schrodinger::rdkit_extensions::isMonomeric(mol)) {
        throw std::invalid_argument("FASTA conversions with atomistic "
                                    "mols are currently unsupported");
    }

    static constexpr std::string_view sru_sgroup{"SRU"};

    std::vector<std::string> errors;
    // monomer repetitions are stored as SRU sgroups
    const auto sgroups = ::RDKit::getSubstanceGroups(mol);
    if (std::any_of(sgroups.begin(), sgroups.end(), [](const auto& sgroup) {
            return get_property<std::string>(&sgroup, SGROUP_TYPE) ==
                   sru_sgroup;
        })) {
        errors.push_back(
            "FASTA conversions with HELMV2.0 monomer repetitions like A'4' and "
            "A'2-4' are currently unsupported.");
    }

    // BLOB and CHEM polymers are unsupported
    auto polymer_ids = get_polymer_ids(mol);
    if (std::any_of(polymer_ids.begin(), polymer_ids.end(),
                    [](const auto& polymer_id) {
                        return !is_biopolymer(polymer_id);
                    })) {
        errors.push_back("FASTA conversions with HELMV2.0 BLOB and "
                         "CHEM polymers are currently unsupported.");
    }

    if (!errors.empty()) {
        throw std::invalid_argument(fmt::format("{}", fmt::join(errors, "\n")));
    }
}

[[nodiscard]] static std::string
get_peptide_fasta(const RDKit::ROMol& mol,
                  const schrodinger::rdkit_extensions::Chain& chain)
{
    using namespace fmt::literals;

    std::string peptide_fasta;
    std::transform(chain.atoms.begin(), chain.atoms.end(),
                   std::back_inserter(peptide_fasta), [&](auto monomer_idx) {
                       auto monomer_id =
                           get_monomer_id(mol.getAtomWithIdx(monomer_idx));
                       auto fasta_monomer =
                           fasta::get_helm_to_fasta_amino_acid(monomer_id);
                       if (!fasta_monomer) {
                           throw std::invalid_argument(fmt::format(
                               "Unsupported monomer: '{}'", monomer_id));
                       }
                       return *fasta_monomer;
                   });

    return fmt::format(">{annotation}\n{monomers}",
                       "annotation"_a = chain.annotation,
                       "monomers"_a = peptide_fasta);
}

[[nodiscard]] static std::string
get_nucleotide_fasta(const RDKit::ROMol& mol,
                     const schrodinger::rdkit_extensions::Chain& chain)
{
    using namespace fmt::literals;

    check_nucleotide_structure(mol, chain);

    std::string nucleotide_fasta;
    for (size_t i = 0; i < chain.atoms.size(); i += 3) {
        auto helm_nucleotide =
            get_monomer_id(mol.getAtomWithIdx(chain.atoms[i + 1]));

        auto fasta_monomer =
            fasta::get_helm_to_fasta_nucleotide(helm_nucleotide);
        if (!fasta_monomer) {
            throw std::invalid_argument(
                fmt::format("Unsupported monomer: '{}'", helm_nucleotide));
        }

        nucleotide_fasta.push_back(*fasta_monomer);
    }

    return fmt::format(">{annotation}\n{monomers}",
                       "annotation"_a = chain.annotation,
                       "monomers"_a = nucleotide_fasta);
}

static void
check_nucleotide_structure(const RDKit::ROMol& mol,
                           const schrodinger::rdkit_extensions::Chain& chain)
{
    if (chain.atoms.size() % 3 != 0) {
        throw std::invalid_argument(
            "FASTA conversions with nucleotides expect all monomers to form "
            "part of a nucleotide unit i.e. R(?)P or [dR](?)P");
    }

    auto reference_sugar = get_property<std::string>(
        mol.getAtomWithIdx(chain.atoms[0]), ATOM_LABEL);
    if (reference_sugar != "R" && reference_sugar != "dR") {
        throw std::invalid_argument(
            "FASTA conversions with nucleotide sugars that aren't 'R' or 'dR' "
            "are currently unsupported.");
    }

    for (size_t i = 0; i < chain.atoms.size(); i += 3) {
        auto sugar = get_property<std::string>(
            mol.getAtomWithIdx(chain.atoms[i]), ATOM_LABEL);
        if (sugar != reference_sugar) {
            throw std::invalid_argument(fmt::format(
                "FASTA conversions with nucleotides must have a uniform sugar "
                "composition along a nucleotide chain."));
        }

        auto phosphate = get_property<std::string>(
            mol.getAtomWithIdx(chain.atoms[i + 2]), ATOM_LABEL);
        if (phosphate != "P") {
            throw std::invalid_argument(
                "FASTA conversions with nucleotides expect the nucleotide "
                "chain to be composed of subunits that contain the sugar, base "
                "and phosphate.");
        }

        // we can do this because there are two bonds in the nucleotide subunit
        // and another bond to the next nucleotide subunit. bonds with indices
        // 0, 3, 6, ... should be the first bond in the nucleotide subunit.
        auto bond_to_base = mol.getBondWithIdx(chain.bonds[i]);
        if (get_property<std::string>(bond_to_base, LINKAGE) !=
            BRANCH_LINKAGE) {
            throw std::invalid_argument(
                "FASTA conversions with a nucleotide subunit that doesn't have "
                "nucleotide base aren't supported");
        }
    }
}
