#define BOOST_TEST_MODULE test_monomer_database

#include <memory>
#include <string_view>

#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"
#include "schrodinger/test/scopedenv.h"
#include "test_common.h"

using namespace schrodinger::rdkit_extensions;

BOOST_DATA_TEST_CASE(Check_Core_Monomers_Are_Canonical,
                     boost::unit_test::data::make({0, 1}), ff_state)
{
    // If this test fails, then we need to update
    // monomer_db_default_monomers.h with new
    // canonical SMILES/CORE_SMILES.
    // (this should only happen if RDKit's canonicalization
    // code changes).

    // We run this with both RDKit's stereo algorithms to check that
    // the monomer DB is agnostic to the flag (we force it enabled).
    // Using ScopedFeatureFlag here won't work, since when this is
    // executed, the toplevel.py script has already set the env var,
    // so just changing the feature flag here won't have any effect. We
    // have to go for the env var directly.
    schrodinger::test::ScopedEnvVar rdk_stereo_algo_env_var(
        "RDK_USE_LEGACY_STEREO_PERCEPTION",
        std::to_string(static_cast<int>(ff_state)));

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

BOOST_AUTO_TEST_CASE(TestNameEscaping)
{
    constexpr std::string_view json_for_12ddR = R"JSON(
        [{
        "symbol": "12ddR",
        "polymer_type": "RNA",
        "natural_analog": "R",
        "smiles": "[H:1]OC[C@H]1OCC[C@@H]1O[H:2]",
         "core_smiles": "OC[C@H]1OCC[C@@H]1O",
         "name": "1',2'-Di-Deoxy-Ribose",
         "monomer_type": "Backbone",
         "author": "Pistoia Alliance",
         "pdbcode": "UNK "
         }]
    )JSON";
    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(json_for_12ddR);
}

BOOST_AUTO_TEST_CASE(TestSmilesEscaping)
{
    constexpr std::string_view json_for_boc_ser_bzl = R"JSON(
        [{
        "symbol": "Boc_Ser(Bzl)",
        "polymer_type": "PEPTIDE",
        "natural_analog": "S",
        "smiles": "CC(C)(C)OC(=O)N[C@@H](COCc1ccccc1)C(=O)[OH:2] |atomProp:7.pdbName. N  :8.pdbName. CA :9.pdbName. CB :10.pdbName. OG :11.pdbName. C' :12.pdbName. C1':13.pdbName. C2':14.pdbName. C3':15.pdbName. C4':16.pdbName. C5':17.pdbName. C6':18.pdbName. C  :19.pdbName. O  |",
         "core_smiles": "CC(C)(C)OC(=O)N[C@H](C=O)COCc1ccccc1",
         "name": "N-alpha-t.-Boc-O-benzyl-L-serine",
         "monomer_type": "Backbone",
         "author": "ChEMBL",
         "pdbcode": "UNK "
         }]
    )JSON";
    auto& monomer_db = MonomerDatabase::instance();
    monomer_db.loadMonomersFromJson(json_for_boc_ser_bzl);
    monomer_db.canonicalizeSmilesFields();
}

BOOST_AUTO_TEST_CASE(TestGetNonNaturalAnalogs)
{
    auto& mdb = MonomerDatabase::instance();
    mdb.resetMonomerDefinitions();

    auto analogs = mdb.getMonomersByNaturalAnalog(ChainType::PEPTIDE);

    // The core DB should contain non-natural peptide analogs
    BOOST_REQUIRE(!analogs.empty());

    // D-Alanine ("dA") should be listed under natural analog "A"
    auto it = analogs.find("A");
    BOOST_REQUIRE(it != analogs.end());
    BOOST_REQUIRE(!it->second.empty());

    bool found_dA = false;
    for (const auto& info : it->second) {
        BOOST_REQUIRE(info.symbol.has_value());
        BOOST_REQUIRE(info.natural_analog.has_value());
        // Every entry must have SYMBOL != NATURAL_ANALOG
        BOOST_CHECK_NE(*info.symbol, *info.natural_analog);
        if (*info.symbol == "dA") {
            found_dA = true;
            BOOST_CHECK_EQUAL(*info.natural_analog, "A");
        }
    }
    BOOST_CHECK(found_dA);

    // Verify that standard amino acids (SYMBOL == NATURAL_ANALOG) are excluded
    for (const auto& [natural_analog, entries] : analogs) {
        for (const auto& info : entries) {
            BOOST_CHECK_NE(info.symbol.value_or(""), natural_analog);
        }
    }
}
