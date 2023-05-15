#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_helm_to_rdkit

#include <boost/test/data/test_case.hpp>
#include <sstream>
#include <string>
#include <vector>

#include <GraphMol/AtomIterators.h>
#include <GraphMol/BondIterators.h>
#include <GraphMol/MonomerInfo.h>

#include "schrodinger/rdkit_extensions/helm/constants.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"
#include "schrodinger/test/checkexceptionmsg.h" // TEST_CHECK_EXCEPTION_MSG_SUBSTR

static auto get_atom_label = [](const auto& atom) {
    return atom->template getProp<std::string>(ATOM_LABEL_PROP_NAME);
};

/*
 * During testing we need an easy way to define what monomers belong to a
 * polymer in a way that can be easily enumerated. This helper class provides an
 * easy way to get the atom and bond information of a helm string without
 * parsing it.
 */
class HELMPolymer
{
  public:
    HELMPolymer(const std::string& polymer_id,
                const std::vector<std::string>& monomers) :
        m_monomers(monomers),
        m_polymer_id(polymer_id)
    {
    }

    std::vector<std::string>::const_iterator begin() const
    {
        return m_monomers.begin();
    }

    std::vector<std::string>::const_iterator end() const
    {
        return m_monomers.end();
    }

    std::string id() const
    {
        return m_polymer_id;
    }

    std::string toString() const
    {
        static auto construct_monomer_id = [](const auto& monomer_id) {
            if (monomer_id.size() == 1 || monomer_id == "BEAD") {
                return monomer_id;
            }
            std::stringstream ss;
            ss << "[" << monomer_id << "]";
            return ss.str();
        };

        std::stringstream ss;
        auto monomer = begin();
        ss << id() << "{" << construct_monomer_id(*monomer);
        while (++monomer != end()) {
            ss << "." << construct_monomer_id(*monomer);
        }
        ss << "}";
        return ss.str();
    }

    unsigned int size() const
    {
        return m_monomers.size();
    }

  private:
    std::vector<std::string> m_monomers;
    std::string m_polymer_id;
};

/*
 * A helper class to house information about the polymers in a way that can be
 * easily queried for testing. This allows us to determine what the output graph
 * of the converted ROMol should look like without having to parse the HELM
 * string.
 */
class HELMInfo
{
  public:
    HELMInfo(const std::vector<HELMPolymer>& polymers) : m_polymers(polymers){};

    std::vector<HELMPolymer>::const_iterator begin() const
    {
        return m_polymers.begin();
    }

    std::vector<HELMPolymer>::const_iterator end() const
    {
        return m_polymers.end();
    }

    std::string toString() const
    {
        std::stringstream ss;
        auto polymer = begin();
        ss << (*polymer).toString();
        while (++polymer != end()) {
            ss << "|" << (*polymer).toString();
        }
        ss << "$$$$V2.0";
        return ss.str();
    }
    const std::vector<HELMPolymer>& getPolymers() const
    {
        return m_polymers;
    }

  private:
    std::vector<HELMPolymer> m_polymers;
};

// Declaring this is required for BOOST_TEST to do introspection
namespace std
{
std::ostream& operator<<(std::ostream& os, const HELMPolymer& polymer)
{
    os << polymer.toString();
    return os;
}

std::ostream& operator<<(std::ostream& os, const HELMInfo& helm)
{
    os << helm.toString();
    return os;
}
}; // namespace std

static auto check_residue_properties = [](const auto& atom,
                                          const auto& chain_id,
                                          const auto& residue_name,
                                          const auto& residue_number) {
    const auto residue_info =
        static_cast<const RDKit::AtomPDBResidueInfo*>(atom->getMonomerInfo());
    BOOST_TEST(residue_info->getResidueName() == residue_name);
    BOOST_TEST(residue_info->getResidueNumber() == residue_number);
    BOOST_TEST(residue_info->getChainId() == chain_id);
};

static auto check_converted_bonds = [](const auto& helm_info, const auto& mol) {
    unsigned int num_bonds = 0;
    for (const auto& polymer : helm_info.getPolymers()) {
        num_bonds += (polymer.size() - 1);
    }
    BOOST_REQUIRE(num_bonds == mol->getNumBonds());

    auto bond = mol->beginBonds();
    for (const auto& polymer : helm_info.getPolymers()) {
        auto atom1 = polymer.begin(), atom2 = polymer.begin();
        while (++atom2 != polymer.end()) {
            BOOST_TEST(*atom1 == get_atom_label((*bond)->getBeginAtom()));
            BOOST_TEST(*atom2 == get_atom_label((*bond)->getEndAtom()));
            ++atom1;
            ++bond;
        }
    }
};

static auto check_converted_atoms = [](const auto& helm_info, const auto& mol) {
    unsigned int num_atoms = 0;
    for (const auto& polymer : helm_info.getPolymers()) {
        num_atoms += polymer.size();
    }
    BOOST_REQUIRE(num_atoms == mol->getNumAtoms());

    auto atom = mol->beginAtoms();
    int chain_id = 1;
    for (const auto& polymer : helm_info.getPolymers()) {
        unsigned int residue_number = 1;
        for (const auto& residue_name : polymer) {
            check_residue_properties(*(atom++), std::to_string(chain_id),
                                     residue_name, residue_number++);
        }
        ++chain_id;
    }
};

static auto check_converted_properties = [](const auto& mol) {
    BOOST_TEST(mol->template getProp<bool>(HELM_MODEL_PROP_NAME));
};

static auto check_helm_conversion = [](const auto& helm_info) {
    const auto& mol = helm_to_rdkit(helm_info.toString());
    check_converted_properties(mol);
    check_converted_atoms(helm_info, mol);
    check_converted_bonds(helm_info, mol);
};

namespace bdata = boost::unit_test::data;

// This test is intended to show the supported HELMV2.0 features. See
// subsequent tests for validation
BOOST_DATA_TEST_CASE(
    TestDemo,
    bdata::make(std::vector<std::string>{
        "PEPTIDE1{A}$$$$V2.0",     // peptides
        "CHEM1{A}$$$$V2.0",        // unknown types
        "RNA1{A}$$$$V2.0",         // nucleotides
        "BLOB1{BEAD}$$$$V2.0",     // beads
        "PEPTIDE1{[dH]}$$$$V2.0",  // multi-character residue names
        "PEPTIDE1{A.[dH]}$$$$V2.0" // multiple monomers
    }),
    input_helm)
{
    helm_to_rdkit(input_helm);
}

BOOST_DATA_TEST_CASE(
    TestTODOs,
    bdata::make(std::vector<std::string>{
        "CHEM1{A}|CHEM2{A}$$$$V2.0", // multiple polymers
        "CHEM1{A}$TOKEN$$$V2.0",     // connections
        "CHEM1{A}$$TOKEN$$V2.0",     // polymer groups
        "CHEM1{A}$$$TOKEN$V2.0",     // extended annotations
        "RNA1{R(A)P}$$$$V2.0",       // branch
        "RNA1{A'5'}$$$$V2.0",        // repetitions
        "RNA1{*}$$$$V2.0",           // unknown
        "PEPTIDE1{G.[[*:1]N[C@@H](C=O)C([*:2])=O].C}$$$$V2.0" // inline smiles
    }) ^ bdata::make(std::vector<std::string>{
             ("Parsing failed around position 9:\n"
              "CHEM1{A}|CHEM2{A}$$$$V2.0\n"
              "________^"),
             ("Parsing failed around position 10:\n"
              "CHEM1{A}$TOKEN$$$V2.0\n"
              "_________^"),
             ("Parsing failed around position 11:\n"
              "CHEM1{A}$$TOKEN$$V2.0\n"
              "__________^"),
             ("Parsing failed around position 12:\n"
              "CHEM1{A}$$$TOKEN$V2.0\n"
              "___________^"),
             ("Parsing failed around position 7:\n"
              "RNA1{R(A)P}$$$$V2.0\n"
              "______^"),
             ("Parsing failed around position 7:\n"
              "RNA1{A'5'}$$$$V2.0\n"
              "______^"),
             ("Parsing failed around position 6:\n"
              "RNA1{*}$$$$V2.0\n"
              "_____^"),
             ("Parsing failed around position 13:\n"
              "PEPTIDE1{G.[[*:1]N[C@@H](C=O)C([*:2])=O].C}$$$$V2.0\n"
              "____________^")}),
    input_helm, expected_err_msg)
{
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(helm_to_rdkit(input_helm),
                                    std::invalid_argument, expected_err_msg);
}
BOOST_DATA_TEST_CASE(TestConversionOfInvalidLinearMonomers,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1$$$$V2.0", "UNSUPPORTED1{A.H.D}$$$$V2.0",
                         "CHEM1{K}}$$$$V2.0", "CHEM1{A.H.D}$$$V2.0",
                         "RNA1{A.dH.D}$$$V2.0", "RNA1{A.[dH.D}$$$V2.0",
                         "PEPTIDE1{A.[[[[[[[[dQ].P}$$$$V2.0",
                         "BLOB1{A.H.D}$$$$V2.0", "PEPTIDE1{[dA\n]}$$$$V2.0"

                     }) ^ bdata::make(std::vector<std::string>{
                              ("Parsing failed around position 9:\n"
                               "PEPTIDE1$$$$V2.0\n"
                               "________^"),
                              ("Parsing failed around position 1:\n"
                               "UNSUPPORTED1{A.H.D}$$$$V2.0\n"
                               "^"),
                              ("Parsing failed around position 9:\n"
                               "CHEM1{K}}$$$$V2.0\n"
                               "________^"),
                              ("Parsing failed around position 8:\n"
                               "CHEM1{A.H.D}$$$V2.0\n"
                               "_______^"),
                              ("Parsing failed around position 8:\n"
                               "RNA1{A.dH.D}$$$V2.0\n"
                               "_______^"),
                              ("Parsing failed around position 11:\n"
                               "RNA1{A.[dH.D}$$$V2.0\n"
                               "__________^"),
                              ("Parsing failed around position 13:\n"
                               "PEPTIDE1{A.[[[[[[[[dQ].P}$$$$V2.0\n"
                               "____________^"),
                              ("Parsing failed around position 7:\n"
                               "BLOB1{A.H.D}$$$$V2.0\n"
                               "______^"),
                              ("Parsing failed around position 13:\n"
                               "PEPTIDE1{[dA\n]}$$$$V2.0\n"
                               "____________^")}),
                     input_helm, expected_err_msg)
{
    TEST_CHECK_EXCEPTION_MSG_SUBSTR(helm_to_rdkit(input_helm),
                                    std::invalid_argument, expected_err_msg);
}

BOOST_DATA_TEST_CASE(TestConversionOfLinearMonomers,
                     bdata::make(std::vector<HELMInfo>{
                         {{{"PEPTIDE1", {"A"}}}},
                         {{{"PEPTIDE1", {"A", "H", "D"}}}},
                         {{{"RNA1", {"R", "A"}}}},
                         {{{"RNA1", {"R", "A", "P"}}}},
                         {{{"CHEM1", {"K"}}}}}),
                     helm_info)
{
    check_helm_conversion(helm_info);
}

BOOST_DATA_TEST_CASE(TestConversionOfMulticharacterMonomers,
                     bdata::make(std::vector<HELMInfo>{
                         {{{"PEPTIDE1", {"dA"}}}},
                         {{{"PEPTIDE1", {"A", "dA"}}}},
                         {{{"PEPTIDE1", {"dG", "P", "dF"}}}},
                         {{{"RNA1", {"dG", "P"}}}},
                         {{{"BLOB1", {"BEAD"}}}}}),
                     helm_info)
{
    check_helm_conversion(helm_info);
}
