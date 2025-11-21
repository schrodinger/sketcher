
#define BOOST_TEST_MODULE test_atom_properties

#include <memory>

#include <boost/test/data/test_case.hpp>

#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/rdkit/atom_properties.h"

using schrodinger::rdkit_extensions::Format;

BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::EnhancedStereo)
BOOST_TEST_DONT_PRINT_LOG_VALUE(
    std::optional<schrodinger::sketcher::EnhancedStereo>)
BOOST_TEST_DONT_PRINT_LOG_VALUE(std::nullopt_t)
BOOST_TEST_DONT_PRINT_LOG_VALUE(schrodinger::sketcher::QueryType)

namespace schrodinger
{
namespace sketcher
{

namespace bdata = boost::unit_test::data;

/**
 * Make sure that an atom with the specified SMILES string produces the expected
 * atom properties.
 */
void check_read_properties_smiles(const std::string& smiles,
                                  const AtomProperties& exp_props)
{
    auto mol = rdkit_extensions::to_rdkit(smiles, Format::SMILES);
    auto* atom = mol->getAtomWithIdx(0);
    auto props = read_properties_from_atom(atom);
    BOOST_TEST(*props == exp_props);
}

/**
 * Make sure that an atom with the specified SMARTS string produces the expected
 * atom query properties.
 */
void check_read_properties_smarts(const std::string& smarts,
                                  const AtomQueryProperties& exp_props)
{
    auto mol = rdkit_extensions::to_rdkit(smarts, Format::SMARTS);
    auto* atom = mol->getAtomWithIdx(0);
    auto props = read_properties_from_atom(atom);
    BOOST_TEST(*props == exp_props);
}

/**
 * Make sure that read_properties_from_atom produces the expected properties for
 * a given input SMILES or SMARTS string
 */
BOOST_AUTO_TEST_CASE(test_read_properties_from_atom)
{
    auto exp_props = AtomProperties();
    exp_props.element = Element::N;
    check_read_properties_smiles("N", exp_props);

    exp_props.charge = 1;
    exp_props.unpaired_electrons = 4;
    check_read_properties_smiles("[N+]", exp_props);

    exp_props = AtomProperties();
    exp_props.element = Element::P;
    exp_props.isotope = 32;
    exp_props.unpaired_electrons = 3;
    check_read_properties_smiles("[32P]", exp_props);

    auto exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::SPECIFIC_ELEMENT;
    exp_query_props.element = Element::N;
    exp_query_props.aromaticity = QueryAromaticity::ALIPHATIC;
    check_read_properties_smarts("N", exp_query_props);

    exp_query_props.aromaticity = QueryAromaticity::AROMATIC;
    exp_query_props.charge = 2;
    exp_query_props.ring_count_type = QueryCount::EXACTLY;
    exp_query_props.ring_count_exact_val = 3;
    check_read_properties_smarts("[nR3++]", exp_query_props);

    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::SPECIFIC_ELEMENT;
    exp_query_props.element = Element::C;
    exp_query_props.aromaticity = QueryAromaticity::ALIPHATIC;
    exp_query_props.isotope = 12;
    check_read_properties_smarts("[12C]", exp_query_props);

    // test allowed atom list queries
    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::ALLOWED_LIST;
    exp_query_props.allowed_list = {Element::O, Element::C, Element::N};
    exp_query_props.aromaticity = QueryAromaticity::ALIPHATIC;
    check_read_properties_smarts("[O,C,N]", exp_query_props);

    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::ALLOWED_LIST;
    exp_query_props.allowed_list = {Element::C, Element::O, Element::N};
    exp_query_props.aromaticity = QueryAromaticity::AROMATIC;
    check_read_properties_smarts("[c,o,n]", exp_query_props);

    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::ALLOWED_LIST;
    exp_query_props.allowed_list = {Element::N, Element::O, Element::F};
    exp_query_props.aromaticity = QueryAromaticity::ANY;
    check_read_properties_smarts("[#7,#8,#9]", exp_query_props);

    // test a long list of elements to ensure that we're preserving the element
    // order correctly
    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::ALLOWED_LIST;
    exp_query_props.allowed_list = {Element::FR, Element::CS, Element::RB,
                                    Element::K,  Element::NA, Element::LI,
                                    Element::F,  Element::CL, Element::BR,
                                    Element::I,  Element::AT, Element::TS};
    exp_query_props.aromaticity = QueryAromaticity::ANY;
    check_read_properties_smarts(
        "[#87,#55,#37,#19,#11,#3,#9,#17,#35,#53,#85,#117]", exp_query_props);

    // test disallowed atom list queries
    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::NOT_ALLOWED_LIST;
    exp_query_props.allowed_list = {Element::NA, Element::MG, Element::K,
                                    Element::MN};
    exp_query_props.aromaticity = QueryAromaticity::ANY;
    check_read_properties_smarts("[!#11&!#12&!#19&!#25]", exp_query_props);

    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::NOT_ALLOWED_LIST;
    exp_query_props.allowed_list = {Element::NA, Element::MG, Element::K,
                                    Element::MN};
    exp_query_props.aromaticity = QueryAromaticity::ANY;
    check_read_properties_smarts("[!Na&!Mg&!K&!Mn]", exp_query_props);

    // test a query that can't be parsed into the dialog properties, so it
    // should register as a query of type SMARTS
    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::SMARTS;
    exp_query_props.smarts_query = "[#7,+2]";
    check_read_properties_smarts("[#7,+2]", exp_query_props);

    // test an R-group query
    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::RGROUP;
    exp_query_props.r_group = 3;
    check_read_properties_smarts("* |$_R3;$|", exp_query_props);

    // an empty query should be treated as an AH wildcard (i.e. any atom)
    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::WILDCARD;
    exp_query_props.wildcard = AtomQuery::AH;
    check_read_properties_smarts("*", exp_query_props);

    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::SMARTS;
    exp_query_props.smarts_query = "[R3,r2]";
    check_read_properties_smarts("[R3,r2]", exp_query_props);

    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::SPECIFIC_ELEMENT;
    exp_query_props.total_h_type = QueryCount::POSITIVE;
    check_read_properties_smarts("[#6&!H0]", exp_query_props);
}

/**
 * Make sure that the specified atom properties produce an atom with the correct
 * SMILES or SMARTS string
 */
void check_create_atom(const std::shared_ptr<AbstractAtomProperties> props,
                       const std::string& exp_smiles_or_smarts, Format format)
{
    auto [atom, enh_stereo] = create_atom_with_properties(props);
    BOOST_TEST(enh_stereo == std::nullopt);
    auto mol = RDKit::RWMol();
    mol.addAtom(atom.get());
    auto smiles_or_smarts = rdkit_extensions::to_string(mol, format);
    BOOST_TEST(smiles_or_smarts == exp_smiles_or_smarts);
}

/**
 * Make sure that create_atom_with_properties produces the expected SMILES or
 * SMARTS string for a given properties object
 */
BOOST_AUTO_TEST_CASE(test_create_atom_with_properties)
{
    auto props = std::make_shared<AtomProperties>();
    props->element = Element::N;
    props->charge = 1;
    props->unpaired_electrons = 4;
    check_create_atom(props, "[N+]", Format::SMILES);

    props->isotope = 13;
    check_create_atom(props, "[13N+]", Format::SMILES);

    auto query_props = std::make_shared<AtomQueryProperties>();
    query_props->query_type = QueryType::SPECIFIC_ELEMENT;
    query_props->element = Element::O;
    query_props->aromaticity = QueryAromaticity::AROMATIC;
    query_props->charge = 2;
    query_props->ring_count_type = QueryCount::EXACTLY;
    query_props->ring_count_exact_val = 3;
    check_create_atom(query_props, "[#8&+2&a&R3]", Format::SMARTS);

    // add an isotope to the properties
    query_props->isotope = 19;
    check_create_atom(query_props, "[#8&19*&+2&a&R3]", Format::SMARTS);

    // an allowed list query
    query_props = std::make_shared<AtomQueryProperties>();
    query_props->query_type = QueryType::ALLOWED_LIST;
    query_props->allowed_list = {Element::MG, Element::K};
    check_create_atom(query_props, "[#12,#19]", Format::SMARTS);

    // a disallowed list
    query_props->query_type = QueryType::NOT_ALLOWED_LIST;
    check_create_atom(query_props, "[!#12&!#19]", Format::SMARTS);

    // check combining allowed and disallowed with other properties
    query_props->num_connections = 2;
    check_create_atom(query_props, "[!#12&!#19&X2]", Format::SMARTS);

    query_props->query_type = QueryType::ALLOWED_LIST;
    check_create_atom(query_props, "[#12,#19;X2]", Format::SMARTS);

    query_props->total_h_type = QueryCount::POSITIVE;
    check_create_atom(query_props, "[#12,#19;!H0;X2]", Format::SMARTS);

    // test R-group
    query_props = std::make_shared<AtomQueryProperties>();
    query_props->query_type = QueryType::RGROUP;
    query_props->r_group = 4;
    check_create_atom(query_props, "* |$_R4$|", Format::EXTENDED_SMILES);

    // test a wildcard, which will get converted to a list of elements in SMARTS
    query_props = std::make_shared<AtomQueryProperties>();
    query_props->query_type = QueryType::WILDCARD;
    query_props->wildcard = AtomQuery::X;
    check_create_atom(query_props, "[#9,#17,#35,#53,#85]", Format::SMARTS);

    // test converting to molv3000 correctly sets properties
    auto output = R"MDL(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 1 0 1 0 0
M  V30 BEGIN ATOM
M  V30 1 C 0.000000 0.000000 0.000000 0 SUBST=-1 RBCNT=-1
M  V30 END ATOM
M  V30 BEGIN SGROUP
M  V30 1 DAT 0 ATOMS=(1 1) QUERYTYPE=SMARTSQ QUERYOP== FIELDDATA="[#6&X0&x0]"
M  V30 END SGROUP
M  V30 END CTAB
M  END
$$$$
)MDL";
    query_props = std::make_shared<AtomQueryProperties>();
    query_props->num_connections = 0;
    query_props->ring_bond_count_type = QueryCount::EXACTLY;
    query_props->ring_bond_count_exact_val = 0;
    check_create_atom(query_props, output, Format::MDL_MOLV3000);
}

/**
 * Make sure that a SMILES string is unchanged when fed through
 * read_properties_from_atom and then create_atom_with_properties
 */
BOOST_DATA_TEST_CASE(test_atom_properties_round_trip_smiles,
                     bdata::make(std::vector<std::string>{
                         "[N+]",
                         "[Fe+2]",
                         "[12C]",
                     }),
                     input_smiles)
{
    auto mol = rdkit_extensions::to_rdkit(input_smiles, Format::SMILES);
    auto atom = mol->getAtomWithIdx(0);
    auto props = read_properties_from_atom(atom);
    auto [output_atom, enh_stereo] = create_atom_with_properties(props);
    BOOST_TEST(enh_stereo == std::nullopt);
    auto output_mol = RDKit::RWMol();
    output_mol.addAtom(output_atom.get());
    auto output_smiles = rdkit_extensions::to_string(
        output_mol, rdkit_extensions::Format::SMILES);
    BOOST_TEST(input_smiles == output_smiles);

    auto second_props = read_properties_from_atom(output_mol.getAtomWithIdx(0));
    BOOST_TEST(*props == *second_props);
}

/**
 * Make sure that a SMARTS string is unchanged when fed through
 * read_properties_from_atom and then create_atom_with_properties
 */
BOOST_DATA_TEST_CASE(test_atom_properties_round_trip_smarts,
                     bdata::make(std::vector<std::string>{
                         "[#7&+&A]",
                         "[#26&+2]",
                         "[#6&12*]",
                         "[#8&+2&R3]",
                         "[#9&!H0]",
                         "[#9,#17,#35,#53,#85]",
                         "[#7,+2]",
                         "[R3,r2]",
                         "[#6&x3]",
                     }),
                     input_smarts)
{
    auto mol = rdkit_extensions::to_rdkit(input_smarts, Format::SMARTS);
    auto atom = mol->getAtomWithIdx(0);
    auto props = read_properties_from_atom(atom);
    auto [output_atom, enh_stereo] = create_atom_with_properties(props);
    auto output_mol = RDKit::RWMol();
    output_mol.addAtom(output_atom.get());
    // make sure the valences are calculated, otherwise we'll get an error when
    // trying to read the properties below
    prepare_mol(output_mol);
    auto output_smarts = rdkit_extensions::to_string(
        output_mol, rdkit_extensions::Format::SMARTS);
    BOOST_TEST(input_smarts == output_smarts);

    auto second_props = read_properties_from_atom(output_mol.getAtomWithIdx(0));
    BOOST_TEST(*props == *second_props);
}

/**
 * Ensure that hasAdvancedProperties() and hasPropertiesBeyondQueryType() return
 * true when the appropriate properties are set to non-default values.
 * Note: hasPropertiesBeyondQueryType excludes stereo properties (SKETCH-2487).
 */
BOOST_AUTO_TEST_CASE(test_hasAdvancedProperties)
{
    AtomQueryProperties query_props;
    BOOST_TEST(!query_props.hasPropertiesBeyondQueryType());
    BOOST_TEST(!query_props.hasAdvancedProperties());
    query_props.query_type = QueryType::SPECIFIC_ELEMENT;
    query_props.element = Element::O;
    BOOST_TEST(!query_props.hasPropertiesBeyondQueryType());
    BOOST_TEST(!query_props.hasAdvancedProperties());

    // Setting stereo should NOT trigger hasPropertiesBeyondQueryType
    // (SKETCH-2487)
    query_props.enhanced_stereo =
        EnhancedStereo(RDKit::StereoGroupType::STEREO_ABSOLUTE, 0);
    BOOST_TEST(!query_props.hasPropertiesBeyondQueryType());
    BOOST_TEST(!query_props.hasAdvancedProperties());
    query_props.enhanced_stereo = std::nullopt;

    query_props.charge = 2;
    BOOST_TEST(query_props.hasPropertiesBeyondQueryType());
    BOOST_TEST(!query_props.hasAdvancedProperties());
    query_props.charge = 0;
    BOOST_TEST(query_props.hasPropertiesBeyondQueryType());
    BOOST_TEST(!query_props.hasAdvancedProperties());
    query_props.charge = std::nullopt;
    BOOST_TEST(!query_props.hasPropertiesBeyondQueryType());
    BOOST_TEST(!query_props.hasAdvancedProperties());
    query_props.aromaticity = QueryAromaticity::AROMATIC;
    BOOST_TEST(query_props.hasPropertiesBeyondQueryType());
    BOOST_TEST(query_props.hasAdvancedProperties());
    query_props.aromaticity = QueryAromaticity::ANY;
    BOOST_TEST(!query_props.hasPropertiesBeyondQueryType());
    BOOST_TEST(!query_props.hasAdvancedProperties());
}

} // namespace sketcher
} // namespace schrodinger