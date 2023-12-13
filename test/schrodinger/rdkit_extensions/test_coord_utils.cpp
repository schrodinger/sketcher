#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_coord_utils

#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/test/data/test_case.hpp>

#include <rdkit/GraphMol/Depictor/DepictUtils.h>

#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/coord_utils.h"

namespace schrodinger
{
namespace rdkit_extensions
{

std::string GIANT_BENZENE = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 6 6 0 0 0
M  V30 BEGIN ATOM
M  V30 1 C -2.340000 -0.060000 0.000000 0
M  V30 2 C -1.340000 1.672051 0.000000 0
M  V30 3 C 0.660000 1.672051 0.000000 0
M  V30 4 C -1.340000 -1.792051 0.000000 0
M  V30 5 C 0.660000 -1.792051 0.000000 0
M  V30 6 C 1.660000 -0.060000 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 5 6
M  V30 2 2 3 6
M  V30 3 2 4 5
M  V30 4 1 1 4
M  V30 5 2 1 2
M  V30 6 1 2 3
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
)";

std::string TOO_SMALL_MOLECULE = R"(
     RDKit          2D

  0  0  0  0  0  0  0  0  0  0999 V3000
M  V30 BEGIN CTAB
M  V30 COUNTS 28 32 0 0 0
M  V30 BEGIN ATOM
M  V30 1 N -0.180969 -4.827934 0.000000 0
M  V30 2 C 0.684718 -4.327610 0.000000 0
M  V30 3 N 0.684341 -3.327343 0.000000 0
M  V30 4 C 1.549928 -2.827192 0.000000 0
M  V30 5 N 1.549625 -1.827197 0.000000 0
M  V30 6 C 2.416364 -3.326887 0.000000 0
M  V30 7 N 3.282051 -2.826562 0.000000 0
M  V30 8 C 4.148313 -3.326158 0.000000 0
M  V30 9 C 5.014000 -2.825833 0.000000 0
M  V30 10 N 5.880263 -3.325429 0.000000 0
M  V30 11 C 6.706370 -2.761842 0.000000 0
M  V30 12 C 6.483613 -1.787165 0.000000 0
M  V30 13 C 7.216485 -1.106732 0.000000 0
M  V30 14 C 8.172187 -1.401248 0.000000 0
M  V30 15 C 8.394845 -2.376098 0.000000 0
M  V30 16 C 7.661974 -3.056531 0.000000 0
M  V30 17 C 8.027705 -3.987299 0.000000 0
M  V30 18 C 7.527762 -4.853363 0.000000 0
M  V30 19 C 6.539072 -5.002846 0.000000 0
M  V30 20 C 6.316569 -5.977700 0.000000 0
M  V30 21 C 5.361167 -6.272597 0.000000 0
M  V30 22 C 4.628094 -5.592542 0.000000 0
M  V30 23 C 4.850424 -4.617588 0.000000 0
M  V30 24 C 5.805826 -4.322690 0.000000 0
M  V30 25 C 4.148443 -4.326053 0.000000 0
M  V30 26 N 3.282756 -4.826377 0.000000 0
M  V30 27 C 2.416667 -4.326881 0.000000 0
M  V30 28 N 1.550807 -4.827106 0.000000 0
M  V30 END ATOM
M  V30 BEGIN BOND
M  V30 1 1 1 2
M  V30 2 2 2 3
M  V30 3 1 3 4
M  V30 4 1 4 5
M  V30 5 2 4 6
M  V30 6 1 6 7
M  V30 7 2 7 8
M  V30 8 1 8 9
M  V30 9 1 9 10
M  V30 10 1 10 11
M  V30 11 2 11 12
M  V30 12 1 12 13
M  V30 13 2 13 14
M  V30 14 1 14 15
M  V30 15 2 15 16
M  V30 16 1 11 16
M  V30 17 1 16 17
M  V30 18 2 17 18
M  V30 19 1 18 19
M  V30 20 2 19 20
M  V30 21 1 20 21
M  V30 22 2 21 22
M  V30 23 1 22 23
M  V30 24 2 23 24
M  V30 25 1 10 24
M  V30 26 1 19 24
M  V30 27 1 8 25
M  V30 28 2 25 26
M  V30 29 1 26 27
M  V30 30 1 6 27
M  V30 31 2 27 28
M  V30 32 1 2 28
M  V30 END BOND
M  V30 END CTAB
M  END
$$$$
)";

namespace utf = boost::unit_test;
namespace bdata = boost::unit_test::data;

/**
 * Make sure that rescale_bond_length_if_needed can successfully rescale both
 * GIANT_BENZENE and TOO_SMALL_MOLECULE
 */
BOOST_TEST_DECORATOR(*utf::tolerance(0.001))
BOOST_DATA_TEST_CASE(test_rescale_bond_length_if_needed,
                     bdata::make(std::vector<std::string>{GIANT_BENZENE,
                                                          TOO_SMALL_MOLECULE}) ^
                         bdata::make(std::vector<double>{2.0, 1.0}),
                     mol_text, exp_starting_bond_length)
{
    auto mol = to_rdkit(mol_text);
    auto starting_bond_length = get_most_common_bond_length(*mol);
    BOOST_TEST(starting_bond_length == exp_starting_bond_length);
    rescale_bond_length_if_needed(*mol);
    auto corrected_bond_length = get_most_common_bond_length(*mol);
    BOOST_TEST(corrected_bond_length == RDDepict::BOND_LEN);
}

/**
 * Make sure that get_most_common_bond_length and rescale_bond_length_if_needed
 * properly handle molecules without a conformer or without any bonds
 */
BOOST_AUTO_TEST_CASE(test_rescaling_mol_with_no_bonds)
{
    auto mol = to_rdkit("C.C");
    BOOST_TEST(!mol->getNumConformers());
    BOOST_TEST(get_most_common_bond_length(*mol) == -1.0);
    // make sure that rescale_bond_length_if_needed doesn't throw an exception
    rescale_bond_length_if_needed(*mol);
    BOOST_TEST(get_most_common_bond_length(*mol) == -1.0);
    compute2DCoords(*mol);
    BOOST_TEST(mol->getNumConformers() == 1);
    BOOST_TEST(get_most_common_bond_length(*mol) == -1.0);
    // make sure that rescale_bond_length_if_needed doesn't throw an exception
    rescale_bond_length_if_needed(*mol);
}

} // namespace rdkit_extensions
} // namespace schrodinger
