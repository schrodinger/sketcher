#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_properties_and_bonds

#include <memory>

#include <boost/test/data/test_case.hpp>

#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

#include "../test_common.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/sketcher/rdkit/atom_properties.h"

namespace schrodinger
{
namespace sketcher
{

namespace bdata = boost::unit_test::data;

/**
 * Make sure that read_properties_from_atom produces the expected properties for
 * a given input SMILES or SMARTS string
 */
BOOST_AUTO_TEST_CASE(test_read_properties_from_atom)
{
    auto mol = RDKit::SmilesToMol("N");
    auto* atom = mol->getAtomWithIdx(0);
    auto props = read_properties_from_atom(atom);
    auto exp_props = AtomProperties();
    exp_props.element = Element::N;
    BOOST_TEST(*props == exp_props);

    mol = RDKit::SmilesToMol("[N+]");
    atom = mol->getAtomWithIdx(0);
    props = read_properties_from_atom(atom);
    exp_props.charge = 1;
    exp_props.unpaired_electrons = 4;
    BOOST_TEST(*props == exp_props);

    mol = RDKit::SmilesToMol("[32P]");
    atom = mol->getAtomWithIdx(0);
    props = read_properties_from_atom(atom);
    exp_props = AtomProperties();
    exp_props.element = Element::P;
    exp_props.isotope = 32;
    exp_props.unpaired_electrons = 3;
    BOOST_TEST(*props == exp_props);

    mol = RDKit::SmartsToMol("N");
    atom = mol->getAtomWithIdx(0);
    props = read_properties_from_atom(atom);
    auto exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::SPECIFIC_ELEMENT;
    exp_query_props.element = Element::N;
    exp_query_props.aromaticity = QueryAromaticity::ALIPHATIC;
    BOOST_TEST(*props == exp_query_props);

    mol = RDKit::SmartsToMol("[nR3++]");
    atom = mol->getAtomWithIdx(0);
    props = read_properties_from_atom(atom);
    exp_query_props.aromaticity = QueryAromaticity::AROMATIC;
    exp_query_props.charge = 2;
    exp_query_props.ring_count_type = QueryCount::EXACTLY;
    exp_query_props.ring_count_exact_val = 3;
    BOOST_TEST(*props == exp_query_props);

    mol = RDKit::SmartsToMol("[12C]");
    atom = mol->getAtomWithIdx(0);
    props = read_properties_from_atom(atom);
    exp_query_props = AtomQueryProperties();
    exp_query_props.query_type = QueryType::SPECIFIC_ELEMENT;
    exp_query_props.element = Element::C;
    exp_query_props.aromaticity = QueryAromaticity::ALIPHATIC;
    exp_query_props.isotope = 12;
    BOOST_TEST(*props == exp_query_props);
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
    auto atom = create_atom_with_properties(props);
    auto mol = RDKit::RWMol();
    mol.addAtom(atom.get());
    auto smiles =
        rdkit_extensions::to_string(mol, rdkit_extensions::Format::SMILES);
    BOOST_TEST(smiles == "[N+]");

    props->isotope = 13;
    auto atom2 = create_atom_with_properties(props);
    auto mol2 = RDKit::RWMol();
    mol2.addAtom(atom2.get());
    auto smiles2 =
        rdkit_extensions::to_string(mol2, rdkit_extensions::Format::SMILES);
    BOOST_TEST(smiles2 == "[13N+]");

    auto query_props = std::make_shared<AtomQueryProperties>();
    query_props->query_type = QueryType::SPECIFIC_ELEMENT;
    query_props->element = Element::O;
    query_props->aromaticity = QueryAromaticity::AROMATIC;
    query_props->charge = 2;
    query_props->ring_count_type = QueryCount::EXACTLY;
    query_props->ring_count_exact_val = 3;
    auto query_atom = create_atom_with_properties(query_props);
    auto query_mol = RDKit::RWMol();
    query_mol.addAtom(query_atom.get());
    auto smarts = rdkit_extensions::to_string(query_mol,
                                              rdkit_extensions::Format::SMARTS);
    BOOST_TEST(smarts == "[#8&+2&a&R3]");

    query_props->isotope = 19;
    auto query_atom2 = create_atom_with_properties(query_props);
    auto query_mol2 = RDKit::RWMol();
    query_mol2.addAtom(query_atom2.get());
    auto smarts2 = rdkit_extensions::to_string(
        query_mol2, rdkit_extensions::Format::SMARTS);
    BOOST_TEST(smarts2 == "[#8&19*&+2&a&R3]");
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
    auto mol = RDKit::SmilesToMol(input_smiles);
    auto atom = mol->getAtomWithIdx(0);
    auto props = read_properties_from_atom(atom);
    auto output_atom = create_atom_with_properties(props);
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
                         "[#7&+&A]", "[#26&+2]", "[#6&12*]", "[#8&+2&R3]"}),
                     input_smarts)
{
    auto mol = RDKit::SmartsToMol(input_smarts);
    auto atom = mol->getAtomWithIdx(0);
    auto props = read_properties_from_atom(atom);
    auto output_atom = create_atom_with_properties(props);
    auto output_mol = RDKit::RWMol();
    output_mol.addAtom(output_atom.get());
    auto output_smarts = rdkit_extensions::to_string(
        output_mol, rdkit_extensions::Format::SMARTS);
    BOOST_TEST(input_smarts == output_smarts);

    auto second_props = read_properties_from_atom(output_mol.getAtomWithIdx(0));
    BOOST_TEST(*props == *second_props);
}

} // namespace sketcher
} // namespace schrodinger