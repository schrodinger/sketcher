#include "schrodinger/rdkit_extensions/monomer_database.h"

#include <iostream>
#include <memory>
#include <optional>
#include <ranges>
#include <regex>
#include <stdexcept>
#include <string>
#include <string_view>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/functional/hash.hpp>
#include <boost/json.hpp>
#include <boost/noncopyable.hpp>

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <rdkit/GraphMol/RDKitBase.h>
#include <rdkit/GraphMol/SmilesParse/SmilesParse.h>
#include <rdkit/GraphMol/SmilesParse/SmilesWrite.h>

#include "sqlite3.h"

#include "schrodinger/rdkit_extensions/monomer_mol.h" // ChainType
#include "schrodinger/rdkit_extensions/monomer_db_schema.h"
#include "schrodinger/rdkit_extensions/monomer_db_default_monomers.h"

using managed_db_t = std::unique_ptr<sqlite3, decltype(&sqlite3_close_v2)>;
using managed_stmt_t =
    std::unique_ptr<sqlite3_stmt, decltype(&sqlite3_finalize)>;

namespace
{
using schrodinger::rdkit_extensions::MonomerInfo;

constexpr std::string_view monomer_defs_table{"monomer_definitions"};
constexpr std::string_view analog_column{"NATURAL_ANALOG"};
constexpr std::string_view author_column{"AUTHOR"};
constexpr std::string_view core_column{"CORE_SMILES"};
constexpr std::string_view monomer_type_column{"MONOMER_TYPE"};
constexpr std::string_view name_column{"NAME"};
constexpr std::string_view pdb_code_column{"PDBCODE"};
constexpr std::string_view polymer_type_column{"POLYMER_TYPE"};
constexpr std::string_view smiles_column{"SMILES"};
constexpr std::string_view symbol_column{"SYMBOL"};

// Legacy column names that we need to temporarily support for backwards
// compatibility, no underscores between words
constexpr std::string_view legacy_analog_column{"NATURALANALOG"};
constexpr std::string_view legacy_monomer_type_column{"MONOMERTYPE"};
constexpr std::string_view legacy_polymer_type_column{"POLYMERTYPE"};

inline const char* _sqlite3_column_cstring(sqlite3_stmt* stmt, int idx)
{
    auto result = sqlite3_column_text(stmt, idx);
    return reinterpret_cast<const char*>(result);
}

// Execute a SQL statement on the DB. No results are expected or parsed
void execute_sql_on(sqlite3* db, std::string_view sql)
{
    char* errmsg = nullptr;
    if (auto rc = sqlite3_exec(db, sql.data(), nullptr, nullptr, &errmsg);
        rc != SQLITE_OK) {
        if (errmsg != nullptr) {
            std::string excmsg{errmsg};
            sqlite3_free(errmsg);
            throw std::runtime_error(excmsg);
        }
        throw std::runtime_error("Unknown error when running sql statement.");
    }
}

// Creates a new SQLite DB file. If db_file == nullptr, then we create
// an in-memory DB.
managed_db_t create_empty_monomer_db()
{
    constexpr auto db_flags =
        SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE | SQLITE_OPEN_MEMORY;

    sqlite3* _db = nullptr;
    if (auto rc = sqlite3_open_v2(":memory:", &_db, db_flags, nullptr);
        rc != SQLITE_OK) {
        throw std::runtime_error("Could not create an in-memory monomer DB.");
    }

    managed_db_t db(_db, &sqlite3_close_v2);

    try {
        execute_sql_on(db.get(), create_monomer_definitions_table_sql.data());
    } catch (const std::runtime_error& e) {
        throw std::runtime_error(fmt::format(
            "Could not create the monomer DB schema: {}", e.what()));
    }

    return db;
}

// Create an empty db and populate the "core" table with the default monomers
managed_db_t create_default_monomers_db()
{
    auto db = create_empty_monomer_db();
    try {
        // The MSVC compiler limits the maximum literal string size to 16Kb,
        // So we have to split this into chunks
        for (const auto& chunk : create_default_monomers) {
            execute_sql_on(db.get(), chunk);
        }
    } catch (const std::runtime_error& e) {
        throw std::runtime_error(fmt::format(
            "Could not insert the default monomer definitions into the DB: {}",
            e.what()));
    }

    return db;
}

// Convenience functions to convert values to std::string
// and trim whitespace in the process
inline std::string convert(const boost::json::value& value)
{
    auto ret = value_to<std::string>(value);
    boost::trim(ret);
    return ret;
}
inline std::string convert(const char* value)
{
    std::string ret{value};
    boost::trim(ret);
    return ret;
}

template <class T>
inline void assign_monomer_info(MonomerInfo& m, std::string key, const T& value)
{
    boost::to_upper(key);
    if (key == symbol_column) {
        m.symbol = convert(value);
    } else if (key == polymer_type_column ||
               key == legacy_polymer_type_column) {
        m.polymer_type = convert(value);
    } else if (key == analog_column || key == legacy_analog_column) {
        m.natural_analog = convert(value);
    } else if (key == smiles_column) {
        m.smiles = convert(value);
    } else if (key == core_column) {
        m.core_smiles = convert(value);
    } else if (key == name_column) {
        m.name = convert(value);
    } else if (key == monomer_type_column ||
               key == legacy_monomer_type_column) {
        m.monomer_type = convert(value);
    } else if (key == author_column) {
        m.author = convert(value);
    } else if (key == pdb_code_column) {
        m.pdbcode = convert(value);
    } else {
        auto msg = fmt::format("unknown field: '{}'", key);
        throw std::runtime_error(msg);
    }
}

inline std::vector<MonomerInfo> json_to_monomer_vector(std::string_view json)
{
    boost::json::value parsed = boost::json::parse(json);
    return boost::json::value_to<std::vector<MonomerInfo>>(parsed);
}

template <class value_getter_t>
[[nodiscard]] auto get_info_from_database(sqlite3* custom_monomers_db,
                                          sqlite3* core_monomers_db,
                                          const std::string& sql,
                                          value_getter_t return_val_converter)
{
    sqlite3_stmt* stmt = nullptr;
    for (sqlite3* db : {custom_monomers_db, core_monomers_db}) {
        if (db != nullptr &&
            sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL) == SQLITE_OK) {
            managed_stmt_t statement(stmt, &sqlite3_finalize);

            // NOTE: Return the first result you see
            if (auto rc = sqlite3_step(stmt);
                rc != SQLITE_DONE && rc != SQLITE_OK) {
                return return_val_converter(stmt);
            };
        }
    }

    // this is an easier way to return nullopt from this templated function
    return std::invoke_result_t<value_getter_t, sqlite3_stmt*>{};
}

std::vector<std::string> get_table_fields(sqlite3* db)
{
    constexpr const char* describe_tbl_cmd =
        "PRAGMA table_info(monomer_definitions);";

    auto column_names_getter = [](sqlite3_stmt* stmt) {
        std::vector<std::string> column_names;
        for (auto rc = sqlite3_step(stmt); rc == SQLITE_ROW;
             rc = sqlite3_step(stmt)) {
            // the column name is always the second field in each
            // row returned by table_info
            column_names.emplace_back(_sqlite3_column_cstring(stmt, 1));
        }
        return column_names;
    };

    return get_info_from_database(db, nullptr, describe_tbl_cmd,
                                  column_names_getter);
}

void insert_monomer_info(sqlite3* db, MonomerInfo& m)
{
    // SMILES and NAME must me unique, so we require them to be specified.
    // We also require PolymerType and Symbol combinations to be unique,
    // but we don't require them to be non-empty. This restriction will be
    // enforced by the SQLite schema.
    for (auto& field : {m.name, m.smiles}) {
        if (!field.has_value() || field->empty()) {
            throw std::runtime_error(
                "Monomer definitions must include at least "
                "SMILES and NAME, which can't be empty.");
        }
    }

    auto insert_command = fmt::format(
        "REPLACE INTO monomer_definitions ("
        "{9}, {10}, {11}, {12}, {13}, {14}, {15}, {16}, {17}"
        ") VALUES ("
        "'{0}', '{1}', '{2}', '{3}', '{4}', '{5}', '{6}', '{7}', '{8}'"
        ");",
        // SMILES and NAME have no default, see the check above!
        m.symbol.value_or(""), m.polymer_type.value_or(""),
        m.natural_analog.value_or(""), m.smiles.value(),
        m.core_smiles.value_or(""), m.name.value(), m.monomer_type.value_or(""),
        m.author.value_or(""), m.pdbcode.value_or(""),
        //
        symbol_column, polymer_type_column, analog_column, smiles_column,
        core_column, name_column, monomer_type_column, author_column,
        pdb_code_column);
    execute_sql_on(db, insert_command.data());
}

template <class T> std::string
get_monomer_query_sql(std::string monomer_id,
                      schrodinger::rdkit_extensions::ChainType polymer_type,
                      const T& selected)
{
    static constexpr std::string_view template_{
        "SELECT {0} FROM {5} WHERE {1}='{2}' AND {3}='{4}';"};

    auto type_value = schrodinger::rdkit_extensions::toString(polymer_type);
    boost::algorithm::trim(monomer_id);

    return fmt::format(template_, selected, polymer_type_column, type_value,
                       symbol_column, monomer_id, monomer_defs_table);
}

std::pair<std::string, std::string>
canonicalize_monomer_smiles(const char* smiles)
{
    using namespace RDKit::v2::SmilesParse;
    static const SmilesParserParams read_params{.removeHs = false,
                                                .replacements = {}};
    static const RDKit::SmilesWriteParams write_params;
    constexpr auto cx_flags = RDKit::SmilesWrite::CXSmilesFields::CX_ATOM_PROPS;

    auto mol = MolFromSmiles(smiles, read_params);
    if (mol == nullptr) {
        auto msg = fmt::format(
            "Error: canonicalization failed, failed to parse SMILES '{}'.",
            smiles);
        throw std::runtime_error(msg);
    }

    // we want to keep atom properties for the pdb atom names.
    auto new_smiles = RDKit::MolToCXSmiles(*mol, write_params, cx_flags);

    // Remove the transformation below once
    // https://github.com/rdkit/rdkit/issues/8664 is fixed
    static const std::regex map_number_regex{
        R"regex(\d+\.molAtomMapNumber\.\d+)regex"};
    static const std::regex consecutive_separators{R"regex(::)regex"};
    new_smiles = std::regex_replace(new_smiles, map_number_regex, "");
    while (new_smiles.find("::") != std::string::npos) {
        new_smiles =
            std::regex_replace(new_smiles, consecutive_separators, ":");
    }
    constexpr std::string_view empty_atom_props = " |atomProp:|";
    if (auto pos = new_smiles.find(empty_atom_props);
        pos != std::string::npos) {
        new_smiles = new_smiles.substr(0, pos);
    }

    // after generating the SMILES, remove the "leaving atoms" for
    // the CORE_SMILES field
    mol->beginBatchEdit();
    for (auto atom : mol->atoms()) {
        if (atom->hasProp(RDKit::common_properties::molAtomMapNumber)) {
            mol->removeAtom(atom);
        }
    }
    mol->commitBatchEdit();

    auto core_smiles = RDKit::MolToSmiles(*mol);

    return {new_smiles, core_smiles};
}

void canonicalize_db(sqlite3* db)
{
    constexpr std::string_view sql_select =
        "SELECT id, smiles FROM monomer_definitions ORDER BY id ASC;";
    constexpr std::string_view sql_replace_tpl =
        "UPDATE monomer_definitions SET smiles='{1}', core_smiles='{2}' "
        "WHERE id={0};";

    execute_sql_on(db, "BEGIN TRANSACTION;");

    try {
        sqlite3_stmt* stmt = nullptr;
        if (auto rc =
                sqlite3_prepare_v2(db, sql_select.data(), -1, &stmt, NULL);
            rc != SQLITE_OK) {
            auto msg =
                fmt::format("Error preparing SMILES canonicalization: {}",
                            sqlite3_errstr(rc));
            throw std::runtime_error(msg);
        }

        managed_stmt_t statement(stmt, &sqlite3_finalize);

        for (auto rc = sqlite3_step(stmt); rc == SQLITE_ROW;
             rc = sqlite3_step(stmt)) {

            // column numbers used here are dictated by sql_select
            auto id = sqlite3_column_int(stmt, 0);
            auto smiles = _sqlite3_column_cstring(stmt, 1);
            auto [new_smiles, new_core_smiles] =
                canonicalize_monomer_smiles(smiles);
            auto sql_replace =
                fmt::format(sql_replace_tpl, id, new_smiles, new_core_smiles);
            execute_sql_on(db, sql_replace.c_str());
        }

    } catch (const std::runtime_error& e) {
        execute_sql_on(db, "ROLLBACK;");
        throw e;
    }

    execute_sql_on(db, "COMMIT;");
}

void insert_monomers_from_json(sqlite3* db, std::string_view json)
{
    auto monomers = json_to_monomer_vector(json);
    for (const auto& m : monomers) {
        if (!m.areRequiredFieldsPopulated()) {
            auto msg = fmt::format(
                "Error: Cannot store monomer definition with missing "
                "fields: {{ {0}='{9}', {1}='{10}', {2}='{11}', {3}='{12}', "
                "{4}='{13}', {5}='{14}', {6}='{15}', {7}='{16}', "
                "{8}='{17}' }}",
                //
                symbol_column, polymer_type_column, analog_column,
                smiles_column, core_column, name_column, monomer_type_column,
                author_column, pdb_code_column,
                //
                m.symbol.value_or(""), m.polymer_type.value_or(""),
                m.natural_analog.value_or(""), m.smiles.value_or(""),
                m.core_smiles.value_or(""), m.name.value_or(""),
                m.monomer_type.value_or(""), m.author.value_or(""),
                m.pdbcode.value_or(""));

            throw std::runtime_error(msg);
        }
    }

    // parse the json and insert all rows in the same transaction
    execute_sql_on(db, "BEGIN TRANSACTION;");

    try {
        for (auto& m : monomers) {
            // Make sure the SMILES/CORE_SMILES are canonical
            // based on the current RDKit version.
            auto [new_smiles, new_core_smiles] =
                canonicalize_monomer_smiles(m.smiles.value().c_str());
            m.smiles = new_smiles;
            m.core_smiles = new_core_smiles;

            if (!m.pdbcode.has_value()) {
                m.pdbcode = "UNL";
            }

            insert_monomer_info(db, m);
        }
    } catch (const std::runtime_error& e) {
        execute_sql_on(db, "ROLLBACK;");
        throw e;
    }

    execute_sql_on(db, "COMMIT;");
}
} // namespace

namespace schrodinger
{
namespace rdkit_extensions
{
// this needs to live in the same namespace as MonomerInfo
static MonomerInfo tag_invoke(boost::json::value_to_tag<MonomerInfo>,
                              boost::json::value const& jv)
{
    MonomerInfo m;
    for (auto& [key, value] : jv.as_object()) {
        assign_monomer_info(m, key.data(), value);
    }

    return m;
}

std::optional<std::string> getMonomerDbPath()
{
#ifdef __EMSCRIPTEN__
    throw std::logic_error("Monomers are not yet supported in WASM Sketcher");
#else
    auto custom_db_path = getenv(CUSTOM_MONOMER_DB_PATH_ENV_VAR.data());
    if (custom_db_path) {
        return std::make_optional<std::string>(custom_db_path);
    }
    return std::nullopt;
#endif
}

size_t MonomerInfo::getHash() const
{
    // The combination of polymer_type and symbol must be unique,
    // though both of them may be an empty string.
    static boost::hash<std::pair<std::string, std::string>> hasher;
    auto p = std::make_pair(polymer_type.value_or(""), symbol.value_or(""));

    return hasher(p);
}

bool MonomerInfo::areRequiredFieldsPopulated() const
{
    return symbol.has_value() &&         //
           polymer_type.has_value() &&   //
           natural_analog.has_value() && //
           smiles.has_value() &&         //
           name.has_value() &&           //
           monomer_type.has_value() &&   //
           author.has_value();
}

void MonomerDatabase::loadMonomersFromSql(std::string_view sql)
{
    managed_db_t db = create_empty_monomer_db();

    execute_sql_on(db.get(), sql);

    if (auto missing_fields = check_db(db.get()); !missing_fields.empty()) {
        auto msg = fmt::format("Error loading SQL into DB: the following "
                               "required columns are missing: {}.",
                               fmt::join(missing_fields, ", "));
        throw std::runtime_error(msg);
    }

    if (auto custom_db_path = getMonomerDbPath(); custom_db_path.has_value()) {
        dumpToFile(db.get(), *custom_db_path);
    }

    // Make sure the SMILES/CORE_SMILES are canonical
    // based on the current RDKit version.
    canonicalize_db(db.get());

    swap_custom_monomers_db(db.release());
}

void MonomerDatabase::loadMonomersFromJson(std::string_view json)
{

    managed_db_t db = create_empty_monomer_db();

    insert_monomers_from_json(db.get(), json);

    if (auto custom_db_path = getMonomerDbPath(); custom_db_path.has_value()) {
        dumpToFile(db.get(), *custom_db_path);
    }

    swap_custom_monomers_db(db.release());
}

void MonomerDatabase::insertMonomersFromJson(std::string_view json)
{
    if (m_custom_monomers_db == nullptr) {
        loadMonomersFromJson(json);
    } else {
        insert_monomers_from_json(m_custom_monomers_db, json);
    }
}

void MonomerDatabase::loadMonomersFromSQLiteFile(
    boost::filesystem::path db_file)
{
    static constexpr std::string_view err_msg_template{
        "Problem opening database: '{}'"};
    sqlite3* _db = nullptr;
    if (auto rc = sqlite3_open_v2(db_file.string().c_str(), &_db,
                                  SQLITE_OPEN_READWRITE, nullptr);
        rc != SQLITE_OK) {
        throw std::runtime_error(
            fmt::format(err_msg_template, db_file.string()));
    }
    managed_db_t db(_db, &sqlite3_close_v2);

    if (auto missing_fields = check_db(db.get()); !missing_fields.empty()) {
        auto msg =
            fmt::format("Error loading DB file '{}': the following required "
                        "columns are missing: {}.",
                        db_file.string(), fmt::join(missing_fields, ", "));
        throw std::runtime_error(msg);
    }

    // Make sure the SMILES/CORE_SMILES are canonical
    // based on the current RDKit version.
    canonicalize_db(db.get());

    swap_custom_monomers_db(db.release());
}

MonomerDatabase& MonomerDatabase::instance()
{
    static MonomerDatabase monomer_db;
    return monomer_db;
}

MonomerDatabase::MonomerDatabase() :
    m_core_monomers_db{create_default_monomers_db().release()}
{
    if (auto path = getMonomerDbPath();
        path.has_value() && boost::filesystem::exists(*path)) {
        try {
            loadMonomersFromSQLiteFile(*path);
        } catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
        }
    }
}

MonomerDatabase::~MonomerDatabase()
{
    resetMonomerDefinitions();
    if (m_core_monomers_db) {
        sqlite3_close_v2(m_core_monomers_db);
    }
}

void MonomerDatabase::swap_custom_monomers_db(sqlite3* db)
{
    std::swap(db, m_custom_monomers_db);
    if (db != nullptr) {
        sqlite3_close_v2(db);
    }
}

std::vector<std::string> MonomerDatabase::check_db(sqlite3* db) const
{
    auto ref_columns = get_table_fields(m_core_monomers_db);

    auto db_columns = get_table_fields(db);
    std::vector<std::string> missing_fields;
    for (auto&& col : ref_columns) {
        if (std::ranges::find(db_columns, col) == db_columns.end()) {
            missing_fields.push_back(std::move(col));
        }
    }
    if (!missing_fields.empty()) {
        return missing_fields;
    }

    return {};
}

void MonomerDatabase::dumpToFile(sqlite3* db,
                                 boost::filesystem::path db_file) const
{
    constexpr const char* main_db = "main";
    constexpr auto db_flags = SQLITE_OPEN_READWRITE | SQLITE_OPEN_CREATE;

    sqlite3* _db = nullptr;
    if (auto rc =
            sqlite3_open_v2(db_file.string().c_str(), &_db, db_flags, nullptr);
        rc != SQLITE_OK) {
        throw std::runtime_error("Could not create an in-memory monomer DB.");
    }

    managed_db_t file_db(_db, &sqlite3_close_v2);

    // See Example 1 in https://sqlite.org/backup.html
    if (auto p_backup =
            sqlite3_backup_init(file_db.get(), main_db, db, main_db);
        p_backup != nullptr) {
        auto result = sqlite3_backup_step(p_backup, -1);

        // free the backup object
        sqlite3_backup_finish(p_backup);
        if (result == SQLITE_DONE) {
            return;
        }
    }

    auto errmsg = fmt::format("Error writing custom monomers to disk: {}.",
                              sqlite3_errmsg(file_db.get()));
    throw std::runtime_error(errmsg);
}

void MonomerDatabase::resetMonomerDefinitions()
{
    swap_custom_monomers_db(nullptr);
}

std::vector<std::string> MonomerDatabase::getDbFields() const
{
    return get_table_fields(m_core_monomers_db);
}

opt_string_t MonomerDatabase::getMonomerSmiles(const std::string& monomer_id,
                                               ChainType polymer_type) const
{
    auto sql = get_monomer_query_sql(monomer_id, polymer_type, smiles_column);

    auto smiles_getter = [](sqlite3_stmt* stmt) -> opt_string_t {
        return std::make_optional(_sqlite3_column_cstring(stmt, 0));
    };

    return get_info_from_database(m_custom_monomers_db, m_core_monomers_db, sql,
                                  smiles_getter);
}

opt_string_t MonomerDatabase::getNaturalAnalog(const std::string& monomer_id,
                                               ChainType polymer_type) const
{
    auto sql = get_monomer_query_sql(monomer_id, polymer_type, analog_column);

    auto analog_getter = [](sqlite3_stmt* stmt) -> opt_string_t {
        return std::make_optional(_sqlite3_column_cstring(stmt, 0));
    };

    return get_info_from_database(m_custom_monomers_db, m_core_monomers_db, sql,
                                  analog_getter);
}

[[nodiscard]] MonomerDatabase::helm_info_t
MonomerDatabase::getHelmInfo(const std::string& pdb_code) const
{
    auto get_sql_command = [&]() -> std::string {
        static constexpr std::string_view template_{
            "SELECT {0}, {1}, {2} FROM {5} WHERE {3}='{4}';"};

        return fmt::format(template_, symbol_column, polymer_type_column,
                           smiles_column, pdb_code_column, pdb_code,
                           monomer_defs_table);
    };

    auto helm_info_getter = [](sqlite3_stmt* stmt) -> helm_info_t {
        auto symbol = _sqlite3_column_cstring(stmt, 0);
        auto polymer_type = _sqlite3_column_cstring(stmt, 1);
        auto smiles = _sqlite3_column_cstring(stmt, 2);
        helm_info_t::value_type helm_info{symbol, smiles,
                                          toChainType(polymer_type)};
        return std::make_optional(helm_info);
    };

    return get_info_from_database(m_custom_monomers_db, m_core_monomers_db,
                                  get_sql_command(), helm_info_getter);
}

[[nodiscard]] opt_string_t
MonomerDatabase::getPdbCode(const std::string& monomer_id,
                            ChainType polymer_type) const
{
    auto sql = get_monomer_query_sql(monomer_id, polymer_type, pdb_code_column);

    auto pdb_code_getter = [](sqlite3_stmt* stmt) -> opt_string_t {
        return std::make_optional(_sqlite3_column_cstring(stmt, 0));
    };

    return get_info_from_database(m_custom_monomers_db, m_core_monomers_db, sql,
                                  pdb_code_getter);
}

[[nodiscard]] std::string
MonomerDatabase::getMonomerDefinitions(bool include_core) const
{
    static constexpr std::string_view template_{
        "SELECT {0} FROM {1} ORDER BY id ASC;"};

    // We don't want "SELECT *", as that would include the "ID" column,
    // which is meaningless and not part of the MonomerInfo schema.
    auto columns = getDbFields();
    const auto sql =
        fmt::format(template_, fmt::join(columns, ", "), monomer_defs_table);

    boost::json::array monomers;

    sqlite3_stmt* stmt = nullptr;
    for (sqlite3* db : {include_core == true ? m_core_monomers_db : nullptr,
                        m_custom_monomers_db}) {
        if (db != nullptr &&
            sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, NULL) == SQLITE_OK) {
            managed_stmt_t statement(stmt, &sqlite3_finalize);

            for (auto rc = sqlite3_step(stmt); rc == SQLITE_ROW;
                 rc = sqlite3_step(stmt)) {

                boost::json::object obj;
                for (int i = 0; i < sqlite3_column_count(stmt); ++i) {
                    auto key = sqlite3_column_name(stmt, i);
                    obj[key] = _sqlite3_column_cstring(stmt, i);
                }
                monomers.push_back(obj);
            }
        }
    }

    return boost::json::serialize(monomers);
}

MonomerInfo MonomerDatabase::getMonomerInfo(const std::string& monomer_id,
                                            ChainType polymer_type) const
{
    // We don't want "SELECT *", as that would include the "ID" column,
    // which is meaningless and not part of the MonomerInfo schema.
    auto fields = fmt::format("{}", fmt::join(getDbFields(), ", "));
    auto sql = get_monomer_query_sql(monomer_id, polymer_type, fields);

    auto monomer_info_getter =
        [](sqlite3_stmt* stmt) -> std::optional<MonomerInfo> {
        auto m = std::make_optional<MonomerInfo>();
        for (int i = 0; i < sqlite3_column_count(stmt); ++i) {
            auto key = sqlite3_column_name(stmt, i);
            auto value = _sqlite3_column_cstring(stmt, i);
            assign_monomer_info(*m, key, value);
        }
        return m;
    };

    if (auto ret = get_info_from_database(
            m_custom_monomers_db, m_core_monomers_db, sql, monomer_info_getter);
        ret.has_value()) {
        return *ret;
    }

    // if not found, return a MonomerInfo with just the type and symbol
    // set (symbol is trimmed). This combination of fields is what we
    // use to generate hashes, and each pair must be unique in the DB.
    return MonomerInfo{
        .symbol = boost::trim_copy(monomer_id),
        .polymer_type = toString(polymer_type),
    };
}

void MonomerDatabase::canonicalizeSmilesFields(bool include_core)
{
    for (sqlite3* db : {include_core == true ? m_core_monomers_db : nullptr,
                        m_custom_monomers_db}) {
        if (db != nullptr) {
            canonicalize_db(db);
        }
    }
}

[[nodiscard]] std::vector<std::pair<std::string, std::string>>
MonomerDatabase::getAllSMILES() const
{
    static constexpr const char* sql =
        "SELECT SMILES, SYMBOL FROM monomer_definitions WHERE "
        "POLYMER_TYPE='PEPTIDE' ORDER BY id ASC;";

    std::vector<std::pair<std::string, std::string>> ret;

    for (sqlite3* db : {m_core_monomers_db, m_custom_monomers_db}) {
        if (db == nullptr) {
            continue;
        }
        sqlite3_stmt* stmt = nullptr;
        if (sqlite3_prepare_v2(db, sql, -1, &stmt, NULL) == SQLITE_OK) {
            managed_stmt_t statement(stmt, &sqlite3_finalize);

            for (auto rc = sqlite3_step(stmt); rc == SQLITE_ROW;
                 rc = sqlite3_step(stmt)) {

                auto smiles = _sqlite3_column_cstring(stmt, 0);
                auto symbol = _sqlite3_column_cstring(stmt, 1);

                ret.emplace_back(smiles, symbol);
            }
        }
    }

    return ret;
}

} // namespace rdkit_extensions
} // namespace schrodinger
