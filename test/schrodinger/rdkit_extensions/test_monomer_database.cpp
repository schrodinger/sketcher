#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test_monomer_database

#include <memory>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"

using namespace schrodinger::rdkit_extensions;

class MonomerDbFixture
{
  public:
    MonomerDbFixture()
    {
        auto db_path = getMonomerDbPath();
        BOOST_REQUIRE(db_path.has_value());

        mdb.reset(new MonomerDatabase(*db_path));
    }

    std::unique_ptr<MonomerDatabase> mdb;
};

BOOST_FIXTURE_TEST_SUITE(TestMonomerDb, MonomerDbFixture);

BOOST_AUTO_TEST_CASE(TestExistingResidue)
{

    constexpr std::string_view ala_smiles =
        "C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB :1.pdbName. CA "
        ":2.pdbName. N  :3.pdbName. H  :4.pdbName. C  :5.pdbName. O  "
        ":6.pdbName. OXT|";

    auto helm_info = mdb->getHelmInfo("ALA");
    BOOST_REQUIRE(helm_info.has_value());
    BOOST_CHECK_EQUAL(std::get<0>(*helm_info), "A");
    BOOST_CHECK_EQUAL(std::get<1>(*helm_info), ala_smiles);
    BOOST_CHECK(std::get<2>(*helm_info) == ChainType::PEPTIDE);

    auto smiles = mdb->getMonomerSmiles("A", ChainType::PEPTIDE);
    BOOST_REQUIRE(smiles.has_value());
    BOOST_CHECK_EQUAL(*smiles, ala_smiles);

    auto pdb_code = mdb->getPdbCode("A", ChainType::PEPTIDE);
    BOOST_REQUIRE(pdb_code.has_value());
    BOOST_CHECK_EQUAL(*pdb_code, "ALA");
}

BOOST_AUTO_TEST_CASE(TestNonExistingResidue)
{
    auto helm_info = mdb->getHelmInfo("DUMMY");
    BOOST_REQUIRE(helm_info.has_value() == false);

    auto smiles = mdb->getMonomerSmiles("DUMMY", ChainType::PEPTIDE);
    BOOST_REQUIRE(smiles.has_value() == false);

    auto pdb_code = mdb->getPdbCode("DUMMY", ChainType::PEPTIDE);
    BOOST_REQUIRE(pdb_code.has_value() == false);
}

BOOST_AUTO_TEST_CASE(TestGetNaturalAnalog)
{
    // SHARED-11627
    // getNaturalAnalog() does NOT return std::optional, like the rest
    // of the MonomerDatabase methods.

    BOOST_CHECK_EQUAL(mdb->getNaturalAnalog("A", ChainType::PEPTIDE), "A");

    BOOST_CHECK_EQUAL(mdb->getNaturalAnalog("DUMMY", ChainType::PEPTIDE), "X");

    // peptide caps do not have analogs!
    BOOST_CHECK_EQUAL(mdb->getNaturalAnalog("ACE", ChainType::PEPTIDE), "X");
}

BOOST_AUTO_TEST_SUITE_END()
