
#define BOOST_TEST_MODULE test_helm_to_rdkit

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include <fmt/format.h>
#include <sstream>
#include <string>
#include <vector>

#include <rdkit/GraphMol/AtomIterators.h>
#include <rdkit/GraphMol/BondIterators.h>
#include <rdkit/GraphMol/MonomerInfo.h>
#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>
#include <rdkit/GraphMol/SubstanceGroup.h>
#include <rdkit/RDGeneral/Invariant.h>

#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"
#include "schrodinger/rdkit_extensions/helm/to_string.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/test/fixtures.h" // silence_stdlog

using namespace helm;

static auto get_atom_label = [](const auto& atom) {
    if (atom->hasProp(MONOMER_LIST)) {
        return atom->template getProp<std::string>(MONOMER_LIST);
    }
    return atom->template getProp<std::string>(ATOM_LABEL);
};

static auto get_sgroups_from_rwmol = [](const auto& rwmol) {
    return ::RDKit::getSubstanceGroups(rwmol);
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
            } else if (monomer_id.find_first_of(",+") != std::string::npos) {
                return fmt::format("({})", monomer_id);
            } else {
                return fmt::format("[{}]", monomer_id);
            }
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
    for (const auto& polymer : helm_info.getPolymers()) {
        unsigned int residue_number = 1;
        for (const auto& residue_name : polymer) {
            check_residue_properties(*(atom++), polymer.id(), residue_name,
                                     residue_number++);
        }
    }
};

static auto check_converted_properties = [](const auto& mol) {
    BOOST_TEST(mol->hasProp(HELM_MODEL));
};

static void check_helm_conversion(const HELMInfo& helm_info)
{
    constexpr auto do_throw = true;
    const auto& mol = helm_to_rdkit(helm_info.toString(), do_throw);
    check_converted_properties(mol);
    check_converted_atoms(helm_info, mol);
    check_converted_bonds(helm_info, mol);
};

namespace bdata = boost::unit_test::data;

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

BOOST_DATA_TEST_CASE(
    TestConversionOfInlineSmiles,
    bdata::make(std::vector<HELMInfo>{
        {{{"PEPTIDE1", {"[*:1]N[C@@H](C=O)C([*:2])=O", "A"}}}},
        {{{"PEPTIDE1", {"A", "[*]N[C@@H](C=O)C([*])=O |$_R1;;;;;;_R2;$|"}}}},
    }),
    helm_info)
{
    check_helm_conversion(helm_info);
}

BOOST_DATA_TEST_CASE(TestConversionOfMonomerLists,
                     bdata::make(std::vector<HELMInfo>{
                         {{{"PEPTIDE1", {"X,G,L"}}}},
                         {{{"PEPTIDE1", {"A", "P+K+H"}}}},
                     }),
                     helm_info)
{
    // NOTE: Not adding monomer lists with ratios because the ratios are
    // discarded during parsing
    check_helm_conversion(helm_info);
}

BOOST_DATA_TEST_CASE(TestConversionOfMultiplePolymers,
                     bdata::make(std::vector<HELMInfo>{
                         {{{"PEPTIDE1", {"dA"}}, {"PEPTIDE2", {"A", "dA"}}}},
                         {{{"PEPTIDE1", {"dG", "X,H,K", "dF"}},
                           {"RNA1", {"dG", "P"}}}},
                         {{{"RNA1", {"dG", "P"}}, {"BLOB1", {"BEAD"}}}}}),
                     helm_info)
{
    check_helm_conversion(helm_info);
}

BOOST_DATA_TEST_CASE(TestMonomerRepetitions,
                     bdata::make(std::vector<std::string>{
                         "A'3'",
                         "(X,C,L)'3-5'",
                         "(A.C.K.L)'3'",
                         "(A.C.K.L)'3-5'",
                     }) ^ bdata::make(std::vector<int>{3, 3, 6, 6}),
                     repeated_monomers, num_monomers)
{
    const auto mol =
        helm_to_rdkit("PEPTIDE1{" + repeated_monomers + "}$$$$V2.0");

    BOOST_TEST(mol->getNumAtoms() == num_monomers);

    if (repeated_monomers.find("-") != std::string::npos) {
        const auto& sgroups = get_sgroups_from_rwmol(*mol);
        const auto& sru_sgroup = *std::find_if(
            sgroups.begin(), sgroups.end(), [](const auto& sgroup) {
                return sgroup.template getProp<std::string>("TYPE") == "SRU";
            });
        // NOTE: not doing a check for the sru_group because test introspection
        // will want me to create a << operator
        BOOST_TEST(sru_sgroup.getAtoms().size() ==
                   static_cast<int>(mol->getNumAtoms() - 2));
        BOOST_TEST(
            (sru_sgroup.template getProp<std::string>("LABEL") == "3" ||
             sru_sgroup.template getProp<std::string>("LABEL") == "3-5"));
    }
};

BOOST_DATA_TEST_CASE(TestConnections,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1,BLOB1,1:R3-?:?",
                         "PEPTIDE1,BLOB1,K:R3-?:?",
                         "PEPTIDE1,BLOB1,1:R3-?:?",
                     }),
                     connections)
{
    const auto mol =
        helm_to_rdkit("PEPTIDE1{K.C}|BLOB1{BEAD}$" + connections + "$$$V2.0");

    // check that there are inter-polymer bonds
    auto num_intra_polymer_bonds = 0;
    for (const auto& sgroup : get_sgroups_from_rwmol(*mol)) {
        if (sgroup.template getProp<std::string>("TYPE") != "COP") {
            continue;
        }
        num_intra_polymer_bonds += sgroup.getBonds().size();
    }
    BOOST_TEST(mol->getNumBonds() > num_intra_polymer_bonds);

    // check that the inter-polymer bonds have the right linkage
    auto bond = ++mol->beginBonds(); // advance from K-C bond
    for (; bond != mol->endBonds(); ++bond) {
        BOOST_TEST((*bond)->template getProp<std::string>(LINKAGE) == "R3-?");
    }
};

BOOST_DATA_TEST_CASE(TestMonomerInlineAnnotations,
                     bdata::make(std::vector<std::string>{
                         R"(A"Something")",
                         "(X,K,L)\"Something\"",
                         R"(R"Something"(C"Something")P"Something")",
                     }),
                     annotated_monomers)
{
    // NOTE: For repeated monomer sequences, we strip off inline
    // annotations if they apply to a repeated subsequence with multiple
    // monomers. e.g., (A.B.C)'N' will have it's annotations stripped
    const auto mol =
        helm_to_rdkit("PEPTIDE1{" + annotated_monomers + "}$$$$V2.0");
    auto atom = mol->beginAtoms();
    for (; atom != mol->endAtoms(); ++atom) {
        // NOTE: to handle dummy query atoms for terminal repeated units
        if (!(*atom)->hasProp(ATOM_LABEL)) {
            continue;
        }
        BOOST_TEST((*atom)->template getProp<std::string>(ANNOTATION) ==
                   "Something");
    }
};

BOOST_DATA_TEST_CASE(TestConnectionInlineAnnotations,
                     bdata::make(std::vector<std::string>{
                         R"(PEPTIDE1,BLOB1,1:R3-?:?"Something")",
                         R"(PEPTIDE1,BLOB1,K:R3-?:?"Something")",
                     }),
                     annotated_connections)
{
    const auto mol = helm_to_rdkit("PEPTIDE1{K.C}|BLOB1{BEAD}$" +
                                   annotated_connections + "$$$V2.0");
    auto bond = mol->beginBonds();
    BOOST_TEST(!(*bond)->hasProp(ANNOTATION)); // K-C bond in the peptide
    ++bond;
    for (; bond != mol->endBonds(); ++bond) {
        BOOST_TEST((*bond)->template getProp<std::string>(ANNOTATION) ==
                   "Something");
    }
};

BOOST_DATA_TEST_CASE(TestPolymerGroups,
                     bdata::make(std::vector<std::string>{
                         "G1(CHEM1,CHEM2)",              // exclusive list
                         "G1(CHEM1+CHEM2)",              // union
                         "G1(CHEM1+CHEM2)|G2(G1,CHEM2)", // nested
                     }),
                     polymer_groups)
{
    // NOTE: Not checking polymer groups with ratios
    // because we currently strip all ratios off polymer
    // list items e.g.,
    // CHEM1{*}|CHEM2{*}$$G1(CHEM1:0.5,CHEM2:0.1-1.3)$$V2.0
    const auto mol =
        helm_to_rdkit("CHEM1{*}|CHEM2{*}$$" + polymer_groups + "$$V2.0");
    auto polymer_group_sgroup = get_sgroups_from_rwmol(*mol).back();
    BOOST_TEST(polymer_groups ==
               polymer_group_sgroup.template getProp<std::vector<std::string>>(
                   "DATAFIELDS")[0]);
};

BOOST_FIXTURE_TEST_CASE(TestExtendedAnnotations,
                        schrodinger::test::silence_stdlog)
{
    // check invalid json
    std::string annotations = "{some invalid json here}";
    BOOST_CHECK_THROW(
        std::ignore = helm_to_rdkit("PEPTIDE1{L}$$$" + annotations + "$V2.0"),
        std::invalid_argument);

    // this should pass
    annotations = R"({"PEPTIDE1":{"ChainType":"hc"}})";
    const auto mol =
        helm_to_rdkit("PEPTIDE1{L.V.A}$$$" + annotations + "$V2.0");

    auto annotation_sgroup = get_sgroups_from_rwmol(*mol).back();
    BOOST_TEST(annotations ==
               annotation_sgroup.template getProp<std::vector<std::string>>(
                   "DATAFIELDS")[1]);
};

BOOST_DATA_TEST_CASE_F(
    schrodinger::test::silence_stdlog, TestInputsWithSelfBond,
    bdata::make(std::vector<std::string>{
        "RNA1{[dR](C)P.[dR](A)P}|RNA2{[dR](G)P.[dR](T)P}$RNA1,RNA1,1:R3-1:R3$$$"
        "V2.0",
        "RNA1{R(C)P.R(A)P}|RNA2{R(G)P.R(U)P}$RNA1,RNA1,1:R3-1:R3$$$V2.0",
    }),
    input_helm)
{
    BOOST_CHECK_THROW(std::ignore = helm_to_rdkit(input_helm),
                      Invar::Invariant);
}

BOOST_DATA_TEST_CASE(TestConversionOfMonomersWithNonstandardNames,
                     bdata::make(std::vector<std::string>{
                         "(N->O)Leu",
                         "(N->O)Val(3-OH)",
                         "-aze",
                         "1-Nal",
                         "2-pyridylmethyl_Gly",
                         "5-Ava",
                         "Abu(5-Tet)",
                         "Ala(O->S)",
                         "Ala(cPent)",
                         "Aoc(2)",
                         "Arg(Me,Me)",
                         "Asp(OMe)",
                         "Asp(Ph(2-NH2))",
                         "Asp_piperidide",
                         "Bal(3-Me)",
                         "Bal(d3-CF3)",
                         "Bn(4-Cl)_Gly",
                         "Cys(EtO2H)_NH2",
                         "Hph(3,4-diCl)",
                         "Hph(4-CF3,3,5-diF)",
                         "Mono21-",
                         "NH2Bu_Gly",
                         "Ser(Ph(2-Cl))",
                         "Sta(3R,4R)",
                         "d(N->O)Gly(allyl)",
                         "dAsp(pyrrol-1-yl)",
                     }) * bdata::make(std::vector<std::string>{
                              "PEPTIDE1{{[{}]}}$$$$V2.0",
                              "PEPTIDE1{{A.[{}].D}}$$$$V2.0",
                              "PEPTIDE1{{A([{}])D}}$$$$V2.0",
                              "RNA1{{R.[{}]}}$$$$V2.0",
                              "RNA1{{R([{}])}}$$$$V2.0",
                              "RNA1{{R.[{}].P}}$$$$V2.0",
                              "RNA1{{R([{}])P}}$$$$V2.0",
                              "CHEM1{{[{}]}}$$$$V2.0",
                          }),
                     test_monomer_id, helm_template)
{
    const auto test_helm =
        fmt::format(fmt::runtime(helm_template), test_monomer_id);
    const auto mol = helm_to_rdkit(test_helm);
    BOOST_TEST(rdkit_to_helm(*mol) == test_helm);
}

BOOST_DATA_TEST_CASE(TestNeighboringMonomerCustomBonds,
                     bdata::make(std::vector<std::string>{
                         "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R3-2:R3$$$V2.0",
                     }),
                     test_helm)
{
    auto mol = helm_to_rdkit(test_helm);

    // Testing this way since affected bond is not the same
    for (auto& bond : mol->bonds()) {
        if (bond->hasProp(CUSTOM_BOND)) {
            BOOST_TEST(bond->template getProp<std::string>(LINKAGE) !=
                       bond->template getProp<std::string>(CUSTOM_BOND));
        }

        // we should be able to roundtrip input
        BOOST_TEST(rdkit_to_helm(*mol) == test_helm);
    }
}

BOOST_DATA_TEST_CASE_F(
    schrodinger::test::silence_stdlog,
    TestNeighboringMonomerUnsupportedCustomBonds,
    bdata::make(std::vector<std::string>{
        "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:pair-2:pair$$$V2.0",
        "PEPTIDE1{C.C}$PEPTIDE1,PEPTIDE1,1:R2-2:R1$$$V2.0",
        "PEPTIDE1{A.C}|PEPTIDE2{C.P}$PEPTIDE1,PEPTIDE2,1:R2-1:"
        "R1|PEPTIDE1,PEPTIDE2,1:R3-1:R3$$$V2.0",
    }),
    test_helm)
{
    BOOST_CHECK_THROW(std::ignore = helm_to_rdkit(test_helm),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestReadingInlineSmilesSurroundedBySquareBrackets)
{
    const auto mol = helm_to_rdkit("PEPTIDE1{[[H][*:1]]}$$$$V2.0");
    const auto test_atom = mol->getAtomWithIdx(0);
    BOOST_TEST(test_atom->template getProp<std::string>(ATOM_LABEL) ==
               "[H][*:1]");
}

BOOST_DATA_TEST_CASE_F(
    schrodinger::test::silence_stdlog, TestInvalidExpressionsV2,
    bdata::make(std::vector<std::string>{
        "BLOB1{}$$$$V2.0",
        "BLOB1{          }$$$$V2.0",
        "BLOB0{HELLO WORLD}$$$$V2.0",
        "BLOB0D{HELLO WORLD}$$$$V2.0",
        "BLOB1{HELLO WORLD\"my annotation\"}$$$$V2.0",
        "BLOB1{HELLO WORLD}$$$$V3.0",
        "BLOB1{HELLO WORLD}$$$not valid json$V2.0",
        "CHEM1{A.B}$$$$V2.0",
        "CHEM1{[]}$$$$V2.0",
        "CHEM1{[A}$$$$V2.0",
        "CHEM1{[A]}$$$$V2.0",
        "CHEM1{dA]}$$$$V2.0",
        "CHEM1{A]}$$$$V2.0",
        "CHEM1{dA}$$$$V2.0",
        "PEPTIDE1{dA}$$$$V2.0",
        "PEPTIDE1{(dA}$$$$V2.0",
        "PEPTIDE1{dA)}$$$$V2.0",
        "PEPTIDE1{(A,[dA)}$$$$V2.0",
        "PEPTIDE1{(A,dA])}$$$$V2.0",
        "PEPTIDE1{(A,[])}$$$$V2.0",
        "PEPTIDE1{(A,[A])}$$$$V2.0",
        "PEPTIDE1{(A:,[A])}$$$$V2.0",
        "PEPTIDE1{(A,[A]:)}$$$$V2.0",
        "PEPTIDE1{(A,[A]+)}$$$$V2.0",
        "PEPTIDE1{(A,[A]+A)}$$$$V2.0",
        "PEPTIDE1{(+A:9+[aA]+A)}$$$$V2.0",
        "PEPTIDE1{(A:9+[aA]+A+)}$$$$V2.0",
        "PEPTIDE1{(A:9++[aA]+A)}$$$$V2.0",
        "PEPTIDE1{(A:0+[aA]+A)}$$$$V2.0",
        "PEPTIDE1{(A:0.1+[aA]+A:0.0)}$$$$V2.0",
        R"(CHEM1{[Az]'3'}$$$$V2.0)",
        "PEPTIDE1{[dA]''}$$$$V2.0",
        "PEPTIDE1{[dA]'hi'}$$$$V2.0",
        "PEPTIDE1{[dA]'hi-10'}$$$$V2.0",
        "PEPTIDE1{[dA]'1-hi'}$$$$V2.0",
        "PEPTIDE1{[dA]'0'}$$$$V2.0",
        "PEPTIDE1{[dA]'000'}$$$$V2.0",
        "PEPTIDE1{[dA]'0-10'}$$$$V2.0",
        "PEPTIDE1{[dA]'-10'}$$$$V2.0",
        "PEPTIDE1{[dA]'10-'}$$$$V2.0",
        "PEPTIDE1{[dA]'1-3-5'}$$$$V2.0",
        "PEPTIDE1{[dA]'01'}$$$$V2.0",
        "PEPTIDE1{[dA]'01-10'}$$$$V2.0",
        "PEPTIDE1{[dA]'3-05'}$$$$V2.0",
        "PEPTIDE1{[dA]'3-1'}$$$$V2.0",
        "PEPTIDE1{[dA]'3-3'}$$$$V2.0",
        "RNA1{R(AP}$$$$V2.0",
        "RNA1{R(A'3')P}$$$$V2.0",
        "RNA1{R(A.L'3')P}$$$$V2.0",
        "RNA1{R(A.L(C)P)P}$$$$V2.0",
        "RNA1{R(A(C)P)P}$$$$V2.0",
        "PEPTIDE1{A.(G.A.C)''.A}$$$$V2.0",
        R"(PEPTIDE1{A(C(A)L)P'9-10'"hello world".A(C)}$$$$V2.0)",
        R"(PEPTIDE1{A(C.A.L)P'9-10'"hello world".A(C)}$$$$V2.0)",
        R"(PEPTIDE1{A((C.A.L))P'9-10'"hello world".A(C)}$$$$V2.0)",
        "PEPTIDE1{A.(G.A.(L.A.P)'9'.C)'3'.A}$$$$V2.0",
        "PEPTIDE1{A.(G.A.C(P))'3'(C).A}$$$$V2.0",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$V1(RNA1,RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$G01(RNA1,RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$G00(RNA1,RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$GA(RNA1,RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$G1(RNA1)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$G1(RNA1,RNA2,)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$G1(,RNA1,RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$G1(RNA1+RNA2+)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$G1(+RNA1+RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$G1(RNA1+RNA2)|$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}$$|G1(RNA1+RNA2)$$V2.0)",
    }),
    test_helm)
{
    // std::invalid_argument should be thrown
    {
        constexpr auto do_throw = true;
        BOOST_CHECK_THROW(std::ignore = helm_to_rdkit(test_helm, do_throw),
                          std::invalid_argument);
    }

    // log and return nullptr
    {
        constexpr auto do_throw = false;
        BOOST_TEST(helm_to_rdkit(test_helm, do_throw) == nullptr);
    }
}

BOOST_DATA_TEST_CASE(
    TestLinearPolymerSupport,
    bdata::make(std::vector<std::string>{
        R"(BLOB1{BEAD}$$$$V2.0)",
        R"(BLOB1{BEAD}"Animated Polystyrene"$$$$V2.0)",
        R"(BLOB1{BEAD}$$${"key":"value"}$V2.0)",
        R"(BLOB1{Bead}"Aminated Polystyrene"|BLOB2{Something here}$$$$V2.0)",
        R"(PEPTIDE1{L.V.A}$$${"PEPTIDE1":{"ChainType":"hc"}}$V2.0)",
        R"(PEPTIDE1{A.G.A.C.A}$$$$V2.0)",
        "CHEM1{[c1ccccc1]}$$$$V2.0",
        R"(CHEM1{[[*:1]OCCCCC=C]}$$$$V2.0)",
        R"(CHEM1{[[*:2]C(=O)[C@H](C)N([*:1])CC]}$$$$V2.0)",
        R"(PEPTIDE1{A.[C(=O)[C@H](C)N([*:1])C]}$$$$V2.0)",
        R"(PEPTIDE1{A.[[*:2]C(=O)[C@H](C)N([*:1])C]}$$$$V2.0)",
        R"(PEPTIDE1{A.C.D}"Polymer 1"|PEPTIDE2{D.A.C}"Polymer 2"$$$$V2.0)",
        R"(CHEM1{*}$$$$V2.0)",
        R"(CHEM1{[SMCC]}|PEPTIDE1{L.M}$$$$V2.0)",
        R"(CHEM1{[sDBl]}$$$$V2.0)",
        R"(PEPTIDE1{A.*.G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.[*C(=O)[C@H](C)N(*)C |$_R2;;;;;;_R1;;;$].A}$$$$V2.0)",
        R"(PEPTIDE1{A.[*C(=O)[C@H](C)N(*)C |$_R2;;;;;;_R1;;;$]}$$$$V2.0)",
        R"(PEPTIDE1{A.[C(=O)[C@H](C)N([*])C |$;;;;;;;_R1;;;$]}$$$$V2.0)",
        R"(PEPTIDE1{A.[C[C@@H](C=O)N(C)(*) |$;;;;;;_R1$|]}$$$$V2.0)",
        R"(PEPTIDE1{A.[[*]C(=O)[C@H](C)N([*])C |$_R2;;;;;;_R1;;;$].A}$$$$V2.0)",
        R"(PEPTIDE1{A.[[*]C(=O)[C@H](C)N([*])C |$_R2;;;;;;_R1;;;$]}$$$$V2.0)",
        R"(PEPTIDE1{A.[meA].C}$$$$V2.0)",
        R"(PEPTIDE1{L.V.A}|PEPTIDE2{L.V.A}$$$$V2.0)",
        R"(PEPTIDE1{[*C(=O)[C@H](C)N(*)C |$_R2;;;;;;_R1;;;$]}$$$$V2.0)",
        R"(PEPTIDE1{[C(=O)[C@H](C)N(*)C |$_;;;;;;;_R1;;;$]}$$$$V2.0)",
        R"(PEPTIDE1{[C(=O)[C@H](C)N([*])C |$_;;;;;;;_R1;;;$]}$$$$V2.0)",
        R"(PEPTIDE1{[C[C@@H](N*)C(*)=O |$;;;_R1;;_R2;;$|].[O=C(*)[C@@H](C(C)C)N* |$;;_R2;;;;;;_R1;$|].[O=C(*)[C@@H](CC(C)C)N* |$;;_R2;;;;;;;_R1;$|]}$$$$V2.0)",
        R"(PEPTIDE1{*"IL6"}$$$$V2.0)",
        R"(PEPTIDE1{A.A.C"mutation".D.E.E}$$$$V2.0)",
        R"(BLOB1{Bead}"Aminated Polystyrene"|PEPTIDE1{A.G.T}$$$$V2.0)",
        R"(CHEM1{[A6OH]}"Test annotation"$$$$V2.0)",
        R"(PEPTIDE1{A.C.D.D.E}"HC"|PEPTIDE2{G.C.S.S.S.P.K.K.V.K}"LC"$$$$V2.0)",
        "PEPTIDE1{[Phe_3Cl]}$$$$V2.0",
        "PEPTIDE1{[Phe-3Cl]}$$$$V2.0",
        "PEPTIDE1{[D-1Nal]}$$$$V2.0",
        "PEPTIDE1{[D-Phe_4F]}$$$$V2.0",
        "RNA1{R.[Phe_3Cl].P}$$$$V2.0",
        "RNA1{R.[Phe-3Cl].P}$$$$V2.0",
        "RNA1{R.[D-1Nal].P}$$$$V2.0",
        "RNA1{R.[D-Phe_4F].P}$$$$V2.0",
        "CHEM1{[Phe_3Cl]}$$$$V2.0",
        "CHEM1{[Phe-3Cl]}$$$$V2.0",
        "CHEM1{[D-1Nal]}$$$$V2.0",
        "CHEM1{[D-Phe_4F]}$$$$V2.0",
        R"(PEPTIDE1{A.(A+G+C).G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.(A,G,L).G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.(A:1+G:2+C:3).G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.(A:1,G:9).G.C}$$$$V2.0)",
        R"(PEPTIDE1{A.A.(X,*).A.A.C.D.D.E.E}$$$$V2.0)",
        R"(PEPTIDE1{A.A.A.A.(A+G:?).A.A.A.A.(A+G+L).A.A.A.A.A.A.A.A.A.A.A.C.D.D.D.D.D.D.D.D.D.D.D.D.D}$$$$V2.0)",
        R"(PEPTIDE1{A.A.A.A.(A:1.1,G:69.5,W:25.5,[Nal]:3.9).A.A.C.D.D.E.E}$$$$V2.0)",
        R"(PEPTIDE1{[fmoc].D.(V,[dV]).Y.A}$$$$V2.0)",
        R"(PEPTIDE1{A'11'}$$$$V2.0)",
        R"(PEPTIDE1{(A,[dA],C:0.4)'11'}$$$$V2.0)",
        R"(PEPTIDE1{A.G.A.C'5'.A}$$$$V2.0)",
        R"(PEPTIDE1{A.G.A.C.A'5-30'}$$$$V2.0)",
        R"(PEPTIDE1{A.G'70-100'"repeatingMonomer"}$$$$V2.0)",
        R"(RNA1{R(A)P.R(U)P}$$$$V2.0)",
        R"(RNA1{R((A,[dA]))P.R(U)P}$$$$V2.0)",
        R"(RNA1{R((A,[dA]:0.2))P.R(U)P}$$$$V2.0)",
        R"(RNA1{R(A"mutation")P.R(U)P}$$$$V2.0)",
        R"(RNA1{R(A)P.R(C)P.R(G)}|RNA2{R(A)P.R(C)P}$$$$V2.0)",
        R"(RNA1{R(A)P.R(C)P.R(G)P}|RNA2{P.R(C)P.R(G)P.R(T)}$$$$V2.0)",
        R"(PEPTIDE1{A.(G.A.C)'3'.A}$$$$V2.0)",
        R"(PEPTIDE1{A.(G.A.C)'3'"my annotation".A}$$$$V2.0)",
        "PEPTIDE1{A.(G.A.C)'3'.A}$$$$V2.0",
        "RNA1{R(A)P.(R(N)P)'4'.(R(G)P)'3-7'}$$$$V2.0",
        "RNA1{R(A)P.R(C)P.(R(T)P)'2'.R(G)}$$$$V2.0",
        R"(PEPTIDE1{A.(C.A.L)'9-10'"hello world".A(C)}$$$$V2.0)",
        "PEPTIDE1{[Bal]'11'(G)C'9'(C)A'9'(L)}$$$$V2.0",
        "PEPTIDE1{[Bal]'11'(G)(G.C)'9'(C)(A.L)'9'(L)}$$$$V2.0",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}|RNA3{R(A)P}$$G1(RNA1,RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}|RNA3{R(A)P}$$G1(RNA1:3,RNA2:2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}|RNA3{R(A)P}$$G1(RNA1+RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}|RNA3{R(A)P}$$G1(RNA1:3+RNA2:2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}|RNA3{R(A)P}$$G1(RNA1,RNA2)|G2(RNA1,RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}|RNA3{R(A)P}$$G1(RNA1,RNA2)|G2(G1,RNA2,RNA3)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}|RNA3{R(A)P}$$G1(RNA1+RNA2)|G2(RNA1+RNA2)$$V2.0)",
        R"(RNA1{R(C)P}"HC"|RNA2{R(U)P}|RNA3{R(A)P}$$G1(RNA1+RNA2)|G2(G1+RNA2+RNA3)$$V2.0)",
    }),
    test_helm)
{
    // turn on exceptions
    {
        constexpr auto do_throw = true;
        auto mol = helm_to_rdkit(test_helm, do_throw);
        BOOST_TEST(rdkit_to_helm(*mol) == test_helm);
    }

    // turn on logging
    {
        constexpr auto do_throw = false;
        auto mol = helm_to_rdkit(test_helm, do_throw);
        BOOST_TEST(rdkit_to_helm(*mol) == test_helm);
    }
}
