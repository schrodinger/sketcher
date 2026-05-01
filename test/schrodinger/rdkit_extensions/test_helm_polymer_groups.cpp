// -------------------------------------------------------------------------
// Unit tests for HELM polymer group APIs
//
// Copyright Schrodinger LLC, All Rights Reserved.
// -------------------------------------------------------------------------
#define BOOST_TEST_MODULE rdkit_extensions_polymer_groups

#include <boost/test/unit_test.hpp>

#include "schrodinger/rdkit_extensions/polymer_group.h"
#include "schrodinger/rdkit_extensions/helm.h"
#include "schrodinger/rdkit_extensions/helm/to_rdkit.h"
#include "schrodinger/rdkit_extensions/helm/to_string.h"

#include <rdkit/GraphMol/RWMol.h>

using namespace schrodinger::rdkit_extensions;

BOOST_AUTO_TEST_SUITE(TestHelmPolymerGroups)

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupBasic)
{
    // Create a simple molecule with two polymers
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$$$V2.0");

    // Add a union polymer group
    add_polymer_group(*mol, "G1", PolymerGroupType::UNION, {"CHEM1", "CHEM2"});

    auto helm_str = helm::rdkit_to_helm(*mol);
    BOOST_TEST(helm_str == "CHEM1{*}|CHEM2{*}$$G1(CHEM1+CHEM2)$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupExclusiveList)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$$$V2.0");

    // Add an exclusive list polymer group
    add_polymer_group(*mol, "G1", PolymerGroupType::EXCLUSIVE_LIST,
                      {"CHEM1", "CHEM2"});

    auto helm_str = helm::rdkit_to_helm(*mol);
    BOOST_TEST(helm_str == "CHEM1{*}|CHEM2{*}$$G1(CHEM1,CHEM2)$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupWithoutRatios)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$$$V2.0");

    // Add a polymer group without ratios (using overload that doesn't take
    // ratios)
    add_polymer_group(*mol, "G1", PolymerGroupType::UNION, {"CHEM1", "CHEM2"});

    auto helm_str = helm::rdkit_to_helm(*mol);
    BOOST_TEST(helm_str == "CHEM1{*}|CHEM2{*}$$G1(CHEM1+CHEM2)$$V2.0");
    BOOST_TEST(get_polymer_groups_helm_string(*mol) == "G1(CHEM1+CHEM2)");
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupWithRatios)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$$$V2.0");

    // Add a polymer group with ratios (using overload with ratios parameter)
    add_polymer_group(*mol, "G1", PolymerGroupType::UNION, {"CHEM1", "CHEM2"},
                      {"2", "5"});

    auto helm_str = helm::rdkit_to_helm(*mol);
    BOOST_TEST(helm_str == "CHEM1{*}|CHEM2{*}$$G1(CHEM1:2+CHEM2:5)$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupWithMixedRatios)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}|CHEM3{*}$$$$V2.0");

    // Add a polymer group where only some entities have ratios
    add_polymer_group(*mol, "G1", PolymerGroupType::UNION,
                      {"CHEM1", "CHEM2", "CHEM3"}, {"2", "", "3.5"});

    auto helm_str = helm::rdkit_to_helm(*mol);
    BOOST_TEST(helm_str ==
               "CHEM1{*}|CHEM2{*}|CHEM3{*}$$G1(CHEM1:2+CHEM2+CHEM3:3.5)$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestAddNestedPolymerGroups)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}|CHEM3{*}$$$$V2.0");

    // Add first group
    add_polymer_group(*mol, "G1", PolymerGroupType::UNION, {"CHEM1", "CHEM2"});

    // Add nested group that references G1
    add_polymer_group(*mol, "G2", PolymerGroupType::UNION, {"G1", "CHEM3"});

    auto helm_str = helm::rdkit_to_helm(*mol);
    BOOST_TEST(
        helm_str ==
        "CHEM1{*}|CHEM2{*}|CHEM3{*}$$G1(CHEM1+CHEM2)|G2(G1+CHEM3)$$V2.0");
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupMultipleEntities)
{
    auto mol =
        helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}|CHEM3{*}|CHEM4{*}$$$$V2.0");

    // Test with more than 2 entities
    add_polymer_group(*mol, "G1", PolymerGroupType::UNION,
                      {"CHEM1", "CHEM2", "CHEM3", "CHEM4"});

    auto helm_str = helm::rdkit_to_helm(*mol);
    BOOST_TEST(helm_str == "CHEM1{*}|CHEM2{*}|CHEM3{*}|CHEM4{*}$$G1(CHEM1+"
                           "CHEM2+CHEM3+CHEM4)$$V2.0");
}

// Error cases

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupInvalidName)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$$$V2.0");

    // Name doesn't start with 'G'
    BOOST_CHECK_THROW(add_polymer_group(*mol, "INVALID",
                                        PolymerGroupType::UNION,
                                        {"CHEM1", "CHEM2"}),
                      std::invalid_argument);

    // Empty name
    BOOST_CHECK_THROW(add_polymer_group(*mol, "", PolymerGroupType::UNION,
                                        {"CHEM1", "CHEM2"}),
                      std::invalid_argument);

    // Name with non-digits after 'G'
    BOOST_CHECK_THROW(add_polymer_group(*mol, "G1A", PolymerGroupType::UNION,
                                        {"CHEM1", "CHEM2"}),
                      std::invalid_argument);

    // Just 'G' without digits
    BOOST_CHECK_THROW(add_polymer_group(*mol, "G", PolymerGroupType::UNION,
                                        {"CHEM1", "CHEM2"}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupTooFewEntities)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$$$V2.0");

    // Only one entity
    BOOST_CHECK_THROW(
        add_polymer_group(*mol, "G1", PolymerGroupType::UNION, {"CHEM1"}),
        std::invalid_argument);

    // Empty entities
    BOOST_CHECK_THROW(
        add_polymer_group(*mol, "G1", PolymerGroupType::UNION, {}),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupEntityNotFound)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$$$V2.0");

    // Reference to non-existent entity
    BOOST_CHECK_THROW(add_polymer_group(*mol, "G1", PolymerGroupType::UNION,
                                        {"CHEM1", "CHEM999"}),
                      std::invalid_argument);

    // Reference to non-existent group
    BOOST_CHECK_THROW(add_polymer_group(*mol, "G1", PolymerGroupType::UNION,
                                        {"CHEM1", "G99"}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupDuplicateName)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}|CHEM3{*}$$$$V2.0");

    // Add first group
    add_polymer_group(*mol, "G1", PolymerGroupType::UNION, {"CHEM1", "CHEM2"});

    // Try to add another group with the same name
    BOOST_CHECK_THROW(add_polymer_group(*mol, "G1", PolymerGroupType::UNION,
                                        {"CHEM2", "CHEM3"}),
                      std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupConflictWithPolymerID)
{
    // Create molecule where a polymer ID could conflict with group name
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$G1(CHEM1+CHEM2)$$V2.0");

    // This should fail because G1 is already a polymer ID
    BOOST_CHECK_THROW(
        add_polymer_group(*mol, "G1", PolymerGroupType::UNION, {"G1", "CHEM2"}),
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestAddPolymerGroupRatiosSizeMismatch)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}|CHEM3{*}$$$$V2.0");

    // Ratios size doesn't match entities size
    BOOST_CHECK_THROW(
        add_polymer_group(*mol, "G1", PolymerGroupType::UNION,
                          {"CHEM1", "CHEM2", "CHEM3"},
                          {"2", "5"}), // only 2 ratios for 3 entities
        std::invalid_argument);
}

BOOST_AUTO_TEST_CASE(TestGetPolymerGroupsEmpty)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$$$V2.0");

    auto polymer_groups = get_polymer_groups_helm_string(*mol);
    BOOST_TEST(polymer_groups.empty());
}

BOOST_AUTO_TEST_CASE(TestGetPolymerGroupsSingle)
{
    auto mol = helm::helm_to_rdkit("CHEM1{*}|CHEM2{*}$$G1(CHEM1+CHEM2)$$V2.0");

    auto polymer_groups = get_polymer_groups_helm_string(*mol);
    BOOST_TEST(polymer_groups == "G1(CHEM1+CHEM2)");
}

BOOST_AUTO_TEST_CASE(TestGetPolymerGroupsMultiple)
{
    auto mol = helm::helm_to_rdkit(
        "CHEM1{*}|CHEM2{*}|CHEM3{*}$$G1(CHEM1+CHEM2)|G2(CHEM2+CHEM3)$$V2.0");

    auto polymer_groups = get_polymer_groups_helm_string(*mol);
    BOOST_TEST(polymer_groups == "G1(CHEM1+CHEM2)|G2(CHEM2+CHEM3)");
}

BOOST_AUTO_TEST_CASE(TestRoundTripPolymerGroups)
{
    // Test that we can round-trip complex polymer groups
    std::vector<std::string> test_cases = {
        "CHEM1{*}|CHEM2{*}$$G1(CHEM1+CHEM2)$$V2.0",
        "CHEM1{*}|CHEM2{*}$$G1(CHEM1,CHEM2)$$V2.0",
        "CHEM1{*}|CHEM2{*}|CHEM3{*}$$G1(CHEM1+CHEM2)|G2(G1+CHEM3)$$V2.0",
    };

    for (const auto& helm_input : test_cases) {
        auto mol = helm::helm_to_rdkit(helm_input);
        auto helm_output = helm::rdkit_to_helm(*mol);
        BOOST_TEST(helm_output == helm_input);
    }
}

BOOST_AUTO_TEST_SUITE_END()
