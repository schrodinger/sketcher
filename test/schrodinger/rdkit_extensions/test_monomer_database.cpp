#define BOOST_TEST_MODULE test_monomer_database

#include <memory>
#include <string_view>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "test_common.h"

using namespace schrodinger::rdkit_extensions;

BOOST_AUTO_TEST_CASE(Check_Core_Monomers_Are_Canonical)
{
    // If this test fails, then we need to update
    // monomer_db_default_monomers.h with new
    // canonical SMILES/CORE_SMILES.
    // (this should only happen if RDKit's canonicalization
    // code changes)

    auto& mdb = MonomerDatabase::instance();

    // Make sure we don't have custom monomers.
    mdb.resetMonomerDefinitions();

    auto current_core_monomers = mdb.getMonomerDefinitions();

    // now, recanonicalize the core monomers
    constexpr bool include_core_monomers = true;
    mdb.canonicalizeSmilesFields(include_core_monomers);

    auto canonical_core_monomers = mdb.getMonomerDefinitions();

    BOOST_CHECK_EQUAL(current_core_monomers, canonical_core_monomers);
}

BOOST_AUTO_TEST_CASE(TestExistingResidue)
{

    constexpr std::string_view ala_smiles =
        "C[C@H](N[H:1])C(=O)[OH:2] |atomProp:0.pdbName. CB :1.pdbName. CA "
        ":2.pdbName. N  :3.pdbName. H  :4.pdbName. C  :5.pdbName. O  "
        ":6.pdbName. OXT|";

    const auto& mdb = MonomerDatabase::instance();

    auto helm_info = mdb.getHelmInfo("ALA");
    BOOST_REQUIRE(helm_info.has_value());
    BOOST_CHECK_EQUAL(std::get<0>(*helm_info), "A");
    BOOST_CHECK_EQUAL(std::get<1>(*helm_info), ala_smiles);
    BOOST_CHECK(std::get<2>(*helm_info) == ChainType::PEPTIDE);

    auto smiles = mdb.getMonomerSmiles("A", ChainType::PEPTIDE);
    BOOST_REQUIRE(smiles.has_value());
    BOOST_CHECK_EQUAL(*smiles, ala_smiles);

    auto pdb_code = mdb.getPdbCode("A", ChainType::PEPTIDE);
    BOOST_REQUIRE(pdb_code.has_value());
    BOOST_CHECK_EQUAL(*pdb_code, "ALA");
}

BOOST_AUTO_TEST_CASE(TestNonExistingResidue)
{
    const auto& mdb = MonomerDatabase::instance();

    auto helm_info = mdb.getHelmInfo("DUMMY");
    BOOST_REQUIRE(helm_info.has_value() == false);

    auto smiles = mdb.getMonomerSmiles("DUMMY", ChainType::PEPTIDE);
    BOOST_REQUIRE(smiles.has_value() == false);

    auto pdb_code = mdb.getPdbCode("DUMMY", ChainType::PEPTIDE);
    BOOST_REQUIRE(pdb_code.has_value() == false);
}

BOOST_AUTO_TEST_CASE(TestGetNaturalAnalog)
{
    // getNaturalAnalog() now returns std::optional, consistent with other
    // MonomerDatabase methods. This resolves the X/N ambiguity issue.
    const auto& mdb = MonomerDatabase::instance();

    auto analog = mdb.getNaturalAnalog("A", ChainType::PEPTIDE);
    BOOST_REQUIRE(analog.has_value());
    BOOST_CHECK_EQUAL(*analog, "A");

    // Non-existent monomers return nullopt instead of "X"
    auto dummy_analog = mdb.getNaturalAnalog("DUMMY", ChainType::PEPTIDE);
    BOOST_CHECK(!dummy_analog.has_value());

    // peptide caps do not have analogs!
    auto ace_analog = mdb.getNaturalAnalog("ACE", ChainType::PEPTIDE);
    BOOST_CHECK(!ace_analog.has_value());
}

BOOST_AUTO_TEST_CASE(TestGetMonomerDbFromSql)
{
    constexpr std::string_view dummy_sql =
        ("    INSERT INTO monomer_definitions (" //
         "SYMBOL, "                              //
         "POLYMER_TYPE, "                        //
         "NATURAL_ANALOG, "                      //
         "SMILES, "                              //
         "CORE_SMILES, "                         //
         "NAME, "                                //
         "MONOMER_TYPE, "                        //
         "AUTHOR, "                              //
         "PDBCODE"                               //
         "   ) VALUES ("                         //
         "'dummy_symbol', "                      //
         "'PEPTIDE', "                           // This one needs to be real!
         "'dummy_analog1', "                     //
         "'NNNN', "                              // This needs to be parseable
         "'dummy_core_smiles1', "                //
         "'dummy_name1', "                       //
         "'dummy_monomertype1', "                //
         "'dummy_author1', "                     //
         "'dummy_pdbcode1'"                      //
         "   ), ("                               //
         "'C', "                                 // Override the core table!
         "'PEPTIDE', "                           // This one needs to be real!
         "'dummy_analog2', "                     //
         "'CCCC', "                              // This needs to be parseable
         "'dummy_core_smiles2', "                //
         "'dummy_name2', "                       //
         "'dummy_monomertype2', "                //
         "'dummy_author2', "                     //
         "'dummy_pdbcode2'"                      //
         ");");

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromSql(dummy_sql);

    // Default monomers are available from the core_monomers table
    auto a_analog = monomer_db.getNaturalAnalog("A", ChainType::PEPTIDE);
    BOOST_REQUIRE(a_analog.has_value());
    BOOST_CHECK_EQUAL(*a_analog, "A");

    // The "C" definition in the custom monomers table overrides/hides the
    // one from the core monomers table
    auto c_analog = monomer_db.getNaturalAnalog("C", ChainType::PEPTIDE);
    BOOST_REQUIRE(c_analog.has_value());
    BOOST_CHECK_EQUAL(*c_analog, "dummy_analog2");

    // This is the one we just added
    auto dummy_analog =
        monomer_db.getNaturalAnalog("dummy_symbol", ChainType::PEPTIDE);
    BOOST_REQUIRE(dummy_analog.has_value());
    BOOST_CHECK_EQUAL(*dummy_analog, "dummy_analog1");

    BOOST_CHECK_THROW(monomer_db.loadMonomersFromSql("not really sql"),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestGetMonomerDbFromJson)
{
    constexpr std::string_view dummy_json =
        ("[{"
         "\"symbol\": \"dummy_symbol\","
         "\"polymer_type\": \"PEPTIDE\"," // This one needs to be real!
         "\"natural_analog\": \"dummy_analog1\","
         "\"smiles\": \"NNNN\"," // This needs to be parseable
         "\"core_smiles\": \"dummy_core_smiles1\","
         "\"name\": \"dummy_name1\","
         "\"monomer_type\": \"dummy_monomertype1\","
         "\"author\": \"dummy_author1\","
         "\"pdbcode\": \"dummy_pdbcode1\""
         "}, {"
         "\"symbol\": \"C\","             // Override the core table!
         "\"polymer_type\": \"PEPTIDE\"," // This one needs to be real!
         "\"natural_analog\": \"dummy_analog2\","
         "\"smiles\": \"CCCC\"," // This needs to be parseable
         "\"core_smiles\": \"dummy_core_smiles2\","
         "\"name\": \"dummy_name2\","
         "\"monomer_type\": \"dummy_monomertype2\","
         "\"author\": \"dummy_author2\","
         "\"pdbcode\": \"dummy_pdbcode2\""
         "}]");

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(dummy_json);

    // Default monomers are available from the core_monomers table
    auto a_analog = monomer_db.getNaturalAnalog("A", ChainType::PEPTIDE);
    BOOST_REQUIRE(a_analog.has_value());
    BOOST_CHECK_EQUAL(*a_analog, "A");

    // The "C" definition in the custom monomers table overrides/hides the
    // one from the core monomers table
    auto c_analog = monomer_db.getNaturalAnalog("C", ChainType::PEPTIDE);
    BOOST_REQUIRE(c_analog.has_value());
    BOOST_CHECK_EQUAL(*c_analog, "dummy_analog2");

    // This is the one we just added
    auto dummy_analog =
        monomer_db.getNaturalAnalog("dummy_symbol", ChainType::PEPTIDE);
    BOOST_REQUIRE(dummy_analog.has_value());
    BOOST_CHECK_EQUAL(*dummy_analog, "dummy_analog1");

    BOOST_CHECK_THROW(monomer_db.loadMonomersFromJson("not really json"),
                      std::runtime_error);
}

BOOST_AUTO_TEST_CASE(TestUpdateCustomDB)
{
    constexpr std::string_view dummy_json =
        ("[{"
         "\"symbol\": \"dummy_symbol\","
         "\"polymer_type\": \"PEPTIDE\"," // This one needs to be real!
         "\"natural_analog\": \"dummy_analog1\","
         "\"smiles\": \"NNNN\"," // This needs to be parseable
         "\"core_smiles\": \"dummy_core_smiles1\","
         "\"name\": \"dummy_name1\","
         "\"monomer_type\": \"dummy_monomertype1\","
         "\"author\": \"dummy_author1\","
         "\"pdbcode\": \"dummy_pdbcode1\""
         "}]");

    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(dummy_json);

    monomer_db.loadMonomersFromSQLiteFile(
        LocalMonomerDbFixture::test_custom_monomer_db);
}