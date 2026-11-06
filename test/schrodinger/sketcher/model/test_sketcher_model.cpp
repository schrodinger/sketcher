#define BOOST_TEST_MODULE Test_Sketcher
#include <memory>

#include <QGraphicsTextItem>
#include <QSignalSpy>
#include <boost/test/data/test_case.hpp>
#include <boost/test/unit_test.hpp>

#include "../test_common.h"
#include "schrodinger/sketcher/model/sketcher_model.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"
#include "schrodinger/rdkit_extensions/monomer_mol.h"

Q_DECLARE_METATYPE(schrodinger::sketcher::ModelKey);
Q_DECLARE_METATYPE(std::unordered_set<schrodinger::sketcher::ModelKey>);
Q_DECLARE_METATYPE(std::unordered_set<QGraphicsItem*>);

BOOST_GLOBAL_FIXTURE(QApplicationRequiredFixture);

BOOST_TEST_DONT_PRINT_LOG_VALUE(QVariant)

namespace schrodinger
{
namespace sketcher
{

/**
 * Subclass that exposes protected data for testing purposes.
 */
class TestSketcherModel : public SketcherModel
{
  public:
    std::unordered_map<ModelKey, QVariant>& getModelMap()
    {
        return m_model_map;
    }
};

/**
 * Verify that the return values of `get_model_keys()` match the keys in the
 * model dictionary.
 */
BOOST_AUTO_TEST_CASE(model_keys)
{

    TestSketcherModel model;
    auto model_map = model.getModelMap();
    auto all_keys = get_model_keys();
    BOOST_TEST(model_map.size() == all_keys.size());
    for (auto& key : all_keys) {
        bool element_found = true;
        try {
            model_map.at(key);
        } catch (const std::out_of_range&) {
            element_found = false;
        }
        BOOST_TEST(element_found);
    }
}

/**
 * Database monomers have a fixed type and do not change peptide symbol types.
 * Each page stores its own database monomer.
 */
BOOST_AUTO_TEST_CASE(monomer_db_monomer)
{
    SketcherModel model;
    using rdkit_extensions::ChainType;
    using rdkit_extensions::MonomerID;
    BOOST_CHECK(get_monomer_db_monomer_key(MonomerToolType::AMINO_ACID) ==
                ModelKey::AMINO_ACID_DB_MONOMER);
    BOOST_CHECK(get_monomer_db_monomer_key(MonomerToolType::NUCLEIC_ACID) ==
                ModelKey::NUCLEIC_ACID_DB_MONOMER);
    for (auto key :
         {ModelKey::AMINO_ACID_DB_MONOMER, ModelKey::NUCLEIC_ACID_DB_MONOMER}) {
        const auto monomer = model.getValue(key).value<MonomerID>();
        BOOST_TEST(monomer.symbol.empty());
        BOOST_CHECK(monomer.chain_type == ChainType::PEPTIDE);
    }
    BOOST_TEST(model.getMonomerDbMonomer().symbol.empty());

    QSignalSpy changed(&model, &SketcherModel::valuesChanged);
    const MonomerID chem{"chem", ChainType::CHEM};
    const MonomerID rna{"rna", ChainType::RNA};
    model.setValue(ModelKey::AMINO_ACID_DB_MONOMER, chem);
    BOOST_REQUIRE_EQUAL(changed.size(), 1);
    model.setValue(ModelKey::AMINO_ACID_DB_MONOMER, chem);
    BOOST_TEST(changed.size() == 1);
    model.setValue(ModelKey::NUCLEIC_ACID_DB_MONOMER, rna);
    BOOST_TEST(changed.size() == 2);

    // each page stores its own monomer, and the no-argument getter follows
    // MONOMER_TOOL_TYPE
    BOOST_CHECK(model.getMonomerDbMonomer(MonomerToolType::AMINO_ACID) == chem);
    BOOST_CHECK(model.getMonomerDbMonomer(MonomerToolType::NUCLEIC_ACID) ==
                rna);
    BOOST_CHECK(model.getMonomerDbMonomer() == chem);
    model.setValue(ModelKey::MONOMER_TOOL_TYPE, MonomerToolType::NUCLEIC_ACID);
    BOOST_CHECK(model.getMonomerDbMonomer() == rna);

    BOOST_CHECK_THROW(
        model.setValue(ModelKey::AMINO_ACID_DB_MONOMER, QString("chem")),
        std::runtime_error);
    BOOST_CHECK_THROW(model.setValue(ModelKey::AMINO_ACID_SYMBOL, chem),
                      std::runtime_error);
    BOOST_TEST(model.getValueString(ModelKey::AMINO_ACID_SYMBOL) == "A");
}

/**
 * Switch the model to the given page and draw tool, as a side-bar button would
 */
static void use_monomeric_draw_tool(SketcherModel& model, MonomerToolType page,
                                    DrawTool draw_tool)
{
    model.setValues(
        {{ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::MONOMERIC)},
         {ModelKey::MONOMER_TOOL_TYPE, QVariant::fromValue(page)},
         {ModelKey::DRAW_TOOL, QVariant::fromValue(draw_tool)}});
}

/**
 * Each page restores the last monomeric draw tool used on it, independently
 * of the other page and of any atomistic or shared tools used in between
 */
BOOST_AUTO_TEST_CASE(get_draw_tool_to_restore)
{
    using rdkit_extensions::ChainType;
    using rdkit_extensions::MonomerID;
    const auto AA = MonomerToolType::AMINO_ACID;
    const auto NA = MonomerToolType::NUCLEIC_ACID;
    SketcherModel model;
    BOOST_TEST((model.getDrawToolToRestore(AA) == DrawTool::MONOMER));
    BOOST_TEST((model.getDrawToolToRestore(NA) == DrawTool::MONOMER));

    use_monomeric_draw_tool(model, AA, DrawTool::MONOMERIC_CONNECTION);
    use_monomeric_draw_tool(model, NA, DrawTool::MONOMER);
    BOOST_TEST(
        (model.getDrawToolToRestore(AA) == DrawTool::MONOMERIC_CONNECTION));
    BOOST_TEST((model.getDrawToolToRestore(NA) == DrawTool::MONOMER));

    // atomistic and shared tools don't replace a page's monomeric tool
    for (auto draw_tool : {DrawTool::BOND, DrawTool::SELECT, DrawTool::ERASE}) {
        model.setValue(ModelKey::DRAW_TOOL, draw_tool);
    }
    model.setValue(ModelKey::MONOMER_TOOL_TYPE, AA);
    BOOST_TEST(
        (model.getDrawToolToRestore(AA) == DrawTool::MONOMERIC_CONNECTION));
    BOOST_TEST((model.getDrawToolToRestore(NA) == DrawTool::MONOMER));

    // a page restores its database monomer tool only when it has a database
    // monomer
    use_monomeric_draw_tool(model, NA, DrawTool::MONOMER_DB_MONOMER);
    BOOST_TEST((model.getDrawToolToRestore(NA) == DrawTool::MONOMER));
    model.setValue(ModelKey::NUCLEIC_ACID_DB_MONOMER,
                   MonomerID{"chem", ChainType::CHEM});
    BOOST_TEST(
        (model.getDrawToolToRestore(NA) == DrawTool::MONOMER_DB_MONOMER));
    BOOST_TEST(
        (model.getDrawToolToRestore(AA) == DrawTool::MONOMERIC_CONNECTION));

    // reset() forgets every page's tool
    model.reset();
    BOOST_TEST((model.getDrawToolToRestore(AA) == DrawTool::MONOMER));
    BOOST_TEST((model.getDrawToolToRestore(NA) == DrawTool::MONOMER));
}

/**
 * A page restores its custom monomer tool only when the stored custom monomer
 * belongs on that page, since both pages share one stored custom monomer
 */
BOOST_AUTO_TEST_CASE(get_draw_tool_to_restore_custom_monomer)
{
    using rdkit_extensions::ChainType;
    const auto AA = MonomerToolType::AMINO_ACID;
    const auto NA = MonomerToolType::NUCLEIC_ACID;
    auto set_custom_monomer = [](SketcherModel& model, const QString& smiles,
                                 ChainType chain_type) {
        model.setValue(ModelKey::CUSTOM_MONOMER,
                       std::make_pair(smiles, chain_type));
    };
    SketcherModel model;
    // the default custom monomer is empty, so there's nothing to restore
    use_monomeric_draw_tool(model, AA, DrawTool::CUSTOM_MONOMER);
    BOOST_TEST((model.getDrawToolToRestore(AA) == DrawTool::MONOMER));

    set_custom_monomer(model, "C[*:1].[*:2]", ChainType::PEPTIDE);
    BOOST_TEST((model.getDrawToolToRestore(AA) == DrawTool::CUSTOM_MONOMER));

    // a nucleic acid custom monomer replaces the peptide, so the amino acid
    // page can't restore it, but the nucleic acid page can
    use_monomeric_draw_tool(model, NA, DrawTool::CUSTOM_MONOMER);
    set_custom_monomer(model, "N[*:1].[*:2]", ChainType::RNA);
    BOOST_TEST((model.getDrawToolToRestore(AA) == DrawTool::MONOMER));
    BOOST_TEST((model.getDrawToolToRestore(NA) == DrawTool::CUSTOM_MONOMER));

    // a CHEM custom monomer belongs on either page
    set_custom_monomer(model, "O[*:1]", ChainType::CHEM);
    BOOST_TEST((model.getDrawToolToRestore(AA) == DrawTool::CUSTOM_MONOMER));
    BOOST_TEST((model.getDrawToolToRestore(NA) == DrawTool::CUSTOM_MONOMER));
}

/**
 * Verify correct behavior of getter method, setter method, and signals.
 */
BOOST_AUTO_TEST_CASE(get_set_signal)
{

    qRegisterMetaType<ModelKey>("ModelKey");
    qRegisterMetaType<std::unordered_set<ModelKey>>(
        "std::unordered_set<ModelKey>");
    auto var = QVariant::fromValue(ModelKey::ATOM_TOOL);
    BOOST_TEST(var.typeId() >= QMetaType::User);
    BOOST_TEST(static_cast<int>(var.value<ModelKey>()) ==
               static_cast<int>(ModelKey::ATOM_TOOL));

    std::unordered_map<ModelKey, QVariant> new_values = {
        {ModelKey::NEW_STRUCTURES_REPLACE_CONTENT, false},
        {ModelKey::SHOW_LID_LEGEND, true},
        {ModelKey::ALLOW_MULTIPLE_RXNS, true},
        {ModelKey::SHOW_VALENCE_ERRORS, false},
        {ModelKey::COLOR_HETEROATOMS, false},
        {ModelKey::SHOW_STEREO_LABELS, false},
        {ModelKey::USE_IMPLICIT_HYDROGENS, true},
        {ModelKey::SELECTION_TOOL, QVariant::fromValue(SelectionTool::LASSO)},
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::MOVE_ROTATE)},
        {ModelKey::ATOM_TOOL, QVariant::fromValue(AtomTool::QUERY)},
        {ModelKey::BOND_TOOL, QVariant::fromValue(BondTool::DOUBLE)},
        {ModelKey::CHARGE_TOOL, QVariant::fromValue(ChargeTool::INCREASE)},
        {ModelKey::RING_TOOL, QVariant::fromValue(RingTool::CYCLOBUTANE)},
        {ModelKey::ENUMERATION_TOOL,
         QVariant::fromValue(EnumerationTool::EXISTING_RGROUP)},
        {ModelKey::MONOMERIC_CONNECTION_TOOL,
         QVariant::fromValue(MonomericConnectionTool::HBOND)},
        {ModelKey::ELEMENT, QVariant::fromValue(Element::N)},
        {ModelKey::ATOM_QUERY, QVariant::fromValue(AtomQuery::Q)},
        {ModelKey::RGROUP_NUMBER, 2u},
        {ModelKey::RESIDUE_TYPE, QString("foo")},
        {ModelKey::MONOMER_TOOL_TYPE,
         QVariant::fromValue(MonomerToolType::NUCLEIC_ACID)},
        {ModelKey::AMINO_ACID_TOOL, QVariant::fromValue(AminoAcidTool::UNK)},
        {ModelKey::AMINO_ACID_SYMBOL, QString("dA")},
        {ModelKey::NUCLEIC_ACID_TOOL, QVariant::fromValue(NucleicAcidTool::dR)},
        {ModelKey::NUCLEIC_ACID_SYMBOL,
         QVariant::fromValue(NucleicAcidMutation{NucleicAcidTool::C, "5mC"})},
        {ModelKey::RNA_NUCLEOBASE, QVariant::fromValue(StdNucleobase::G)},
        {ModelKey::DNA_NUCLEOBASE, QVariant::fromValue(StdNucleobase::U_OR_T)},
        {ModelKey::CUSTOM_NUCLEOTIDE,
         QVariant::fromValue(std::tuple<std::string, std::string, std::string>(
             "Tho", "I", "PS"))},
        {ModelKey::CUSTOM_MONOMER,
         QVariant::fromValue(
             std::make_pair<QString, rdkit_extensions::ChainType>(
                 "CC", rdkit_extensions::ChainType::CHEM))},
        {ModelKey::INTERFACE_TYPE, InterfaceType::ATOMISTIC},
        {ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::MONOMERIC)},
        {ModelKey::MOLECULE_TYPE,
         QVariant::fromValue(MoleculeType::ATOMISTIC)}};

    TestSketcherModel model;
    QSignalSpy pinged_spy(&model, &SketcherModel::valuePinged);
    QSignalSpy changed_spy(&model, &SketcherModel::valuesChanged);
    std::vector<QSignalSpy*> spies = {&pinged_spy, &changed_spy};
    for (auto& key : get_model_keys()) {
        if (key == ModelKey::RESIDUE_TYPE ||
            key == ModelKey::AMINO_ACID_SYMBOL ||
            key == ModelKey::CUSTOM_MONOMER ||
            key == ModelKey::AMINO_ACID_DB_MONOMER ||
            key == ModelKey::NUCLEIC_ACID_DB_MONOMER ||
            key == ModelKey::NUCLEIC_ACID_SYMBOL ||
            key == ModelKey::CUSTOM_NUCLEOTIDE) {
            // These values are not stored as int-like objects
            continue;
        }
        // Access the current value
        const auto& value = model.getValue(key);
        for (auto spy : spies) {
            BOOST_TEST(spy->empty());
        }
        BOOST_TEST(value.isValid());

        // Assign a new value
        auto new_value = new_values[key];
        model.setValue(key, new_value);

        // Verify that a signal was emitted
        for (auto spy : spies) {
            BOOST_TEST(spy->count() == 1);
        }

        {
            // valuePinged() is emitted with a key and corresponding value
            auto args = pinged_spy.takeLast();
            BOOST_TEST(args.at(0).typeId() >= QMetaType::User);
            auto signal_key = args.at(0).value<ModelKey>();
            BOOST_TEST(static_cast<int>(signal_key) == static_cast<int>(key));
            BOOST_TEST(args.at(1).toInt() == new_value.toInt());
        }

        {
            // valuesChanged() is emitted with a set of changed keys
            auto args = changed_spy.takeLast();
            BOOST_TEST(args.at(0).typeId() >= QMetaType::User);
            auto keys = args.at(0).value<std::unordered_set<ModelKey>>();
            BOOST_TEST(keys.size() == 1);
        }

        // Assign the same value again, and verify that the "pinged" signal is
        // emitted but the "changed" signal is not
        model.setValue(key, new_value);
        BOOST_TEST(pinged_spy.count() == 1);
        auto args = pinged_spy.takeLast();
        BOOST_TEST(args.at(0).typeId() >= QMetaType::User);
        auto signal_key = args.at(0).value<ModelKey>();
        BOOST_TEST(static_cast<int>(signal_key) == static_cast<int>(key));
        BOOST_TEST(args.at(1).toInt() == new_value.toInt());
        BOOST_TEST(changed_spy.empty());
    }
}

/**
 * Verify correct behavior of multiple setter method `setValues()`.
 */
BOOST_AUTO_TEST_CASE(setValues)
{

    qRegisterMetaType<std::unordered_set<ModelKey>>(
        "std::unordered_set<ModelKey>");
    TestSketcherModel model;
    QSignalSpy changed_spy(&model, &SketcherModel::valuesChanged);

    std::unordered_set<ModelKey> keys_to_change = {ModelKey::SELECTION_TOOL,
                                                   ModelKey::DRAW_TOOL};
    std::unordered_set<ModelKey> keys_stay_the_same = {
        ModelKey::ENUMERATION_TOOL, ModelKey::ELEMENT};
    std::unordered_map<ModelKey, QVariant> kv_pairs = {
        {ModelKey::SELECTION_TOOL, QVariant::fromValue(SelectionTool::LASSO)},
        {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::BOND)},
    };
    for (auto key : keys_to_change) {
        BOOST_TEST(model.getValue(key) != kv_pairs[key]);
    }
    for (auto key : keys_stay_the_same) {
        kv_pairs.emplace(key, model.getValue(key));
    }

    BOOST_TEST(kv_pairs.size() ==
               keys_to_change.size() + keys_stay_the_same.size());

    model.setValues(kv_pairs);
    BOOST_TEST(changed_spy.count() == 1);
    auto args = changed_spy.takeLast();
    BOOST_TEST(args.at(0).typeId() >= QMetaType::User);
    auto signal_keys = args.at(0).value<std::unordered_set<ModelKey>>();

    std::function<bool(ModelKey)> is_in_signal_keys =
        [signal_keys](ModelKey key) { return signal_keys.count(key) == 1; };
    BOOST_TEST(std::all_of(keys_to_change.begin(), keys_to_change.end(),
                           is_in_signal_keys));
    BOOST_TEST(std::none_of(keys_stay_the_same.begin(),
                            keys_stay_the_same.end(), is_in_signal_keys));

    for (auto& [key, value] : kv_pairs) {
        BOOST_TEST(model.getValue(key) == value);
    }
}

/**
 * Undo any changes to the monomer database when a test finishes
 */
struct ResetDatabase {
    ~ResetDatabase()
    {
        rdkit_extensions::MonomerDatabase::instance().resetMonomerDefinitions();
    }
};

const rdkit_extensions::MonomerID PEPTIDE_DB_MONOMER{
    "pepX", rdkit_extensions::ChainType::PEPTIDE};
const rdkit_extensions::MonomerID RNA_DB_MONOMER{
    "rnaN", rdkit_extensions::ChainType::RNA};
const rdkit_extensions::MonomerID CHEM_DB_MONOMER{
    "chem", rdkit_extensions::ChainType::CHEM};

/**
 * Load a custom monomer database containing the requested monomers out of
 * PEPTIDE_DB_MONOMER, RNA_DB_MONOMER and CHEM_DB_MONOMER
 */
static void load_db_monomers(bool peptide, bool rna, bool chem)
{
    std::vector<std::string> entries;
    if (peptide) {
        entries.push_back(
            R"({"symbol":"pepX","polymer_type":"PEPTIDE","natural_analog":"X",
                "smiles":"CC","name":"Peptide","monomer_type":"backbone",
                "author":"test"})");
    }
    if (rna) {
        entries.push_back(
            R"({"symbol":"rnaN","polymer_type":"RNA","natural_analog":"N",
                "smiles":"CCC","name":"Base","monomer_type":"branch",
                "author":"test"})");
    }
    if (chem) {
        entries.push_back(
            R"({"symbol":"chem","polymer_type":"CHEM","natural_analog":"",
                "smiles":"CN","name":"Chem","monomer_type":"undefined",
                "author":"test"})");
    }
    std::string json = "[";
    for (std::size_t i = 0; i < entries.size(); ++i) {
        json += (i > 0 ? "," : "") + entries[i];
    }
    json += "]";
    auto result =
        rdkit_extensions::MonomerDatabase::instance().loadMonomersFromJson(
            json);
    BOOST_REQUIRE(result.second.empty());
    BOOST_REQUIRE_EQUAL(result.first.size(), entries.size());
}

/**
 * Store the model values for activating a database monomer on the given page
 * without a canvas selection
 */
static void activate_db_monomer(SketcherModel& model, MonomerToolType page,
                                const rdkit_extensions::MonomerID& monomer)
{
    model.setValues(
        {{ModelKey::DRAW_TOOL,
          QVariant::fromValue(DrawTool::MONOMER_DB_MONOMER)},
         {ModelKey::MONOMER_TOOL_TYPE, QVariant::fromValue(page)},
         {get_monomer_db_monomer_key(page), QVariant::fromValue(monomer)}});
}

/**
 * A database change that keeps the remembered database monomers in the Monomer
 * Database table leaves the model unchanged
 */
BOOST_AUTO_TEST_CASE(monomer_db_repair_keeps_valid_monomers)
{
    ResetDatabase reset_database;
    load_db_monomers(true, true, true);
    SketcherModel model;
    activate_db_monomer(model, MonomerToolType::AMINO_ACID, CHEM_DB_MONOMER);
    activate_db_monomer(model, MonomerToolType::NUCLEIC_ACID, RNA_DB_MONOMER);
    QSignalSpy pinged(&model, &SketcherModel::valuePinged);
    QSignalSpy changed(&model, &SketcherModel::valuesChanged);
    load_db_monomers(false, true, true);
    BOOST_TEST(pinged.empty());
    BOOST_TEST(changed.empty());
    BOOST_CHECK(model.getMonomerDbMonomer(MonomerToolType::AMINO_ACID) ==
                CHEM_DB_MONOMER);
    BOOST_CHECK(model.getMonomerDbMonomer(MonomerToolType::NUCLEIC_ACID) ==
                RNA_DB_MONOMER);
    for (auto page :
         {MonomerToolType::AMINO_ACID, MonomerToolType::NUCLEIC_ACID}) {
        BOOST_TEST(
            (model.getDrawToolToRestore(page) == DrawTool::MONOMER_DB_MONOMER));
    }
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);

    // core monomers with no natural analog are always in the table
    const rdkit_extensions::MonomerID core_monomer{
        "ac", rdkit_extensions::ChainType::PEPTIDE};
    activate_db_monomer(model, MonomerToolType::AMINO_ACID, core_monomer);
    rdkit_extensions::MonomerDatabase::instance().resetMonomerDefinitions();
    BOOST_CHECK(model.getMonomerDbMonomer() == core_monomer);
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
}

/**
 * Removing the active page's database monomer while its tool is active clears
 * the page's last-used database monomer and switches to the page's last
 * standard monomer tool, which activating the database monomer left alone,
 * without pinging
 */
BOOST_DATA_TEST_CASE(monomer_db_repair_active_page,
                     boost::unit_test::data::make({true, false}), amino_acid)
{
    ResetDatabase reset_database;
    load_db_monomers(true, true, true);
    SketcherModel model;
    const auto page = amino_acid ? MonomerToolType::AMINO_ACID
                                 : MonomerToolType::NUCLEIC_ACID;
    // the last standard monomer buttons used: Gly and C
    model.setValues(
        {{ModelKey::AMINO_ACID_TOOL, QVariant::fromValue(AminoAcidTool::GLY)},
         {ModelKey::AMINO_ACID_SYMBOL, QString("G")},
         {ModelKey::NUCLEIC_ACID_TOOL, QVariant::fromValue(NucleicAcidTool::C)},
         {ModelKey::NUCLEIC_ACID_SYMBOL,
          QVariant::fromValue(NucleicAcidMutation{NucleicAcidTool::C, "C"})}});
    activate_db_monomer(model, page,
                        amino_acid ? PEPTIDE_DB_MONOMER : RNA_DB_MONOMER);
    BOOST_TEST(
        (model.getDrawToolToRestore(page) == DrawTool::MONOMER_DB_MONOMER));

    QSignalSpy pinged(&model, &SketcherModel::valuePinged);
    QSignalSpy changed(&model, &SketcherModel::valuesChanged);
    load_db_monomers(!amino_acid, amino_acid, true);
    BOOST_TEST(pinged.empty());
    BOOST_REQUIRE_EQUAL(changed.size(), 1);
    const auto changed_keys =
        changed.front().front().value<std::unordered_set<ModelKey>>();
    BOOST_TEST(changed_keys.contains(get_monomer_db_monomer_key(page)));
    BOOST_TEST(changed_keys.contains(ModelKey::DRAW_TOOL));
    for (auto key :
         {ModelKey::AMINO_ACID_TOOL, ModelKey::AMINO_ACID_SYMBOL,
          ModelKey::NUCLEIC_ACID_TOOL, ModelKey::NUCLEIC_ACID_SYMBOL}) {
        BOOST_TEST(!changed_keys.contains(key));
    }

    BOOST_TEST(model.getMonomerDbMonomer(page).symbol.empty());
    BOOST_TEST(
        (model.getDrawToolToRestore(page) != DrawTool::MONOMER_DB_MONOMER));
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
    BOOST_CHECK(model.getMonomerToolType() == page);
    // both pages keep their last standard monomer buttons
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::GLY);
    BOOST_TEST(model.getValueString(ModelKey::AMINO_ACID_SYMBOL) == "G");
    BOOST_CHECK(model.getNucleicAcidTool() == NucleicAcidTool::C);
    const auto symbol = model.getValue(ModelKey::NUCLEIC_ACID_SYMBOL)
                            .value<NucleicAcidMutation>();
    BOOST_CHECK(symbol.tool == NucleicAcidTool::C);
    BOOST_TEST(symbol.symbol == "C");
}

/**
 * Removing an inactive page's database monomer, or removing a database monomer
 * while atomistic mode is active, clears the page's last-used database monomer
 * without changing the active tool, so the page restores its standard monomer
 * tool instead
 */
BOOST_AUTO_TEST_CASE(monomer_db_repair_inactive_page)
{
    ResetDatabase reset_database;
    load_db_monomers(true, true, true);
    SketcherModel model;
    activate_db_monomer(model, MonomerToolType::NUCLEIC_ACID, RNA_DB_MONOMER);
    activate_db_monomer(model, MonomerToolType::AMINO_ACID, PEPTIDE_DB_MONOMER);

    // remove the inactive nucleic acid page's monomer
    load_db_monomers(true, false, true);
    BOOST_TEST(model.getMonomerDbMonomer(MonomerToolType::NUCLEIC_ACID)
                   .symbol.empty());
    BOOST_CHECK(model.getNucleicAcidTool() == NucleicAcidTool::RNA_NUCLEOTIDE);
    BOOST_TEST((model.getDrawToolToRestore(MonomerToolType::NUCLEIC_ACID) !=
                DrawTool::MONOMER_DB_MONOMER));
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER_DB_MONOMER);
    BOOST_CHECK(model.getMonomerDbMonomer() == PEPTIDE_DB_MONOMER);
    BOOST_TEST((model.getDrawToolToRestore(MonomerToolType::AMINO_ACID) ==
                DrawTool::MONOMER_DB_MONOMER));

    // remove the amino acid page's monomer while atomistic mode is active
    model.setValues(
        {{ModelKey::TOOL_SET, QVariant::fromValue(ToolSet::ATOMISTIC)},
         {ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::BOND)}});
    load_db_monomers(false, false, true);
    BOOST_TEST(model.getMonomerDbMonomer().symbol.empty());
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::ALA);
    BOOST_TEST((model.getDrawToolToRestore(MonomerToolType::AMINO_ACID) !=
                DrawTool::MONOMER_DB_MONOMER));
    BOOST_CHECK(model.getDrawTool() == DrawTool::BOND);
}

/**
 * A page whose last tool was a standard monomer keeps it when its last-used
 * database monomer is removed
 */
BOOST_AUTO_TEST_CASE(monomer_db_repair_keeps_standard_tool)
{
    ResetDatabase reset_database;
    load_db_monomers(true, true, true);
    SketcherModel model;
    // a mutation stores the last-used database monomer without changing the
    // draw tool or the page's standard tool
    model.setValues(
        {{ModelKey::DRAW_TOOL, QVariant::fromValue(DrawTool::SELECT)},
         {ModelKey::AMINO_ACID_TOOL, QVariant::fromValue(AminoAcidTool::CYS)},
         {ModelKey::AMINO_ACID_SYMBOL, QString("C")},
         {ModelKey::AMINO_ACID_DB_MONOMER,
          QVariant::fromValue(PEPTIDE_DB_MONOMER)}});
    load_db_monomers(false, true, true);
    BOOST_TEST(model.getMonomerDbMonomer().symbol.empty());
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::CYS);
    BOOST_TEST(model.getValueString(ModelKey::AMINO_ACID_SYMBOL) == "C");
    BOOST_CHECK(model.getDrawTool() == DrawTool::SELECT);
}

/**
 * Unloading a custom override of a core monomer that has a natural analog
 * removes the monomer from the Monomer Database table, so the page's last-used
 * database monomer is cleared even though the core monomer can still be drawn
 */
BOOST_AUTO_TEST_CASE(monomer_db_repair_unloaded_override)
{
    ResetDatabase reset_database;
    auto& db = rdkit_extensions::MonomerDatabase::instance();
    db.resetMonomerDefinitions();
    auto core_meg =
        db.getMonomerInfo("meG", rdkit_extensions::ChainType::PEPTIDE);
    BOOST_REQUIRE(core_meg.natural_analog == std::string("G"));
    auto result = db.loadMonomersFromJson(
        R"([{"symbol":"meG","polymer_type":"PEPTIDE","natural_analog":"X",
             "smiles":"CC","name":"Custom meG","monomer_type":"backbone",
             "author":"test"}])");
    BOOST_REQUIRE(result.second.empty());
    SketcherModel model;
    const rdkit_extensions::MonomerID meg{"meG",
                                          rdkit_extensions::ChainType::PEPTIDE};
    activate_db_monomer(model, MonomerToolType::AMINO_ACID, meg);
    db.resetMonomerDefinitions();
    BOOST_TEST(model.getMonomerDbMonomer().symbol.empty());
    BOOST_CHECK(model.getAminoAcidTool() == AminoAcidTool::ALA);
    BOOST_CHECK(model.getDrawTool() == DrawTool::MONOMER);
}

} // namespace sketcher
} // namespace schrodinger
