#include "schrodinger/rdkit_extensions/fasta/monomers.h"

#include <algorithm>
#include <boost/assign.hpp>
#include <boost/bimap.hpp>
#include <initializer_list>
#include <unordered_map>

#include <fmt/format.h>

namespace fasta
{

// clang-format off
static const boost::bimap<char, std::string> supported_peptides_bimap =
    boost::assign::list_of<boost::bimap<char, std::string>::relation>
        ('A', "A") // alanine
        ('C', "C") // cysteine
        ('D', "D") // aspartate
        ('E', "E") // glutamate
        ('F', "F") // phenylalanine
        ('G', "G") // glycine
        ('H', "H") // histidine
        ('I', "I") // isoleucine
        ('K', "K") // lysine
        ('L', "L") // leucine
        ('M', "M") // methionine
        ('N', "N") // asparagine
        ('P', "P") // proline
        ('Q', "Q") // glutamine
        ('R', "R") // arginine
        ('S', "S") // serine
        ('T', "T") // threonine
        ('V', "V") // valine
        ('W', "W") // tryptophan
        ('Y', "Y") // tyrosine
        // Non-standard
        ('O', "O") // pyrrolysine
        ('U', "U") // selenocysteine
        // Unknown monomer
        ('X', "X"); // any, meaning one single unknown amino acid

static const boost::bimap<char, std::string> supported_nucleotides_bimap =
    boost::assign::list_of<boost::bimap<char, std::string>::relation>
        ('A', "A")    // adenine
        ('C', "C")    // cytosine
        ('G', "G")    // guanine
        ('T', "T")    // thymine
        ('U', "U")    // uracil
        // Unknown monomer
        ('N', "N"); // any, meaning one single unknown base
// clang-format on

template <class T, class U>
static auto get_monomer_helper(const T& monomer_map, U& query)
{
    auto monomer = monomer_map.find(query);
    return monomer == monomer_map.end() ? std::nullopt
                                        : std::optional(monomer->second);
}

[[nodiscard]] std::optional<std::string>
get_fasta_to_helm_amino_acid(char one_letter_monomer)
{
    return get_monomer_helper(supported_peptides_bimap.left,
                              one_letter_monomer);
};

[[nodiscard]] std::optional<std::string>
get_fasta_to_helm_nucleotide(char one_letter_monomer,
                             std::string_view sugar_helm_monomer)
{
    // NOTE: the helm entries have to separate the phosphate, sugar and
    // nucleotide base subunits
    auto monomer = get_monomer_helper(supported_nucleotides_bimap.left,
                                      one_letter_monomer);
    return monomer == std::nullopt
               ? std::nullopt
               : std::optional<std::string>(
                     fmt::format("{}({})P", sugar_helm_monomer, *monomer));
}

[[nodiscard]] std::optional<char>
get_helm_to_fasta_amino_acid(const std::string& helm_amino_acid)
{
    return get_monomer_helper(supported_peptides_bimap.right, helm_amino_acid);
}

[[nodiscard]] std::optional<char>
get_helm_to_fasta_nucleotide(const std::string& helm_nucleotide)
{
    return get_monomer_helper(supported_nucleotides_bimap.right,
                              helm_nucleotide);
}

} // namespace fasta
