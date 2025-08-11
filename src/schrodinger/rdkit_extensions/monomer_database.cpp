#include "schrodinger/rdkit_extensions/monomer_database.h"

#include <array>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include <boost/algorithm/string/trim.hpp>
#include <boost/noncopyable.hpp>
#include <boost/filesystem.hpp>
#include <fmt/format.h>

#include "schrodinger/rdkit_extensions/monomer_mol.h" // ChainType

#include "sqlite3.h"

using namespace std::string_view_literals;

namespace schrodinger
{
namespace rdkit_extensions
{

std::optional<std::string> getCustomMonomerDbPath()
{
#ifdef __EMSCRIPTEN__
    throw std::logic_error("Monomers are not yet supported in WASM Sketcher");
#else
    auto custom_db_path = getenv(CUSTOM_MONOMER_DB_PATH_ENV_VAR.c_str());
    if (custom_db_path) {
        return std::make_optional<std::string>(custom_db_path);
    }
    return std::nullopt;
#endif
}

std::optional<std::string> getMonomerDbPath()
{
#ifdef __EMSCRIPTEN__
    throw std::logic_error("Monomers are not yet supported in WASM Sketcher");
#else
    auto custom_db_path = getCustomMonomerDbPath();
    if (custom_db_path.has_value() &&
        boost::filesystem::exists(*custom_db_path)) {
        return custom_db_path;
    }
    auto db_path = getenv(DEFAULT_MONOMER_DB_PATH_ENV_VAR.c_str());
    if (db_path) {
        return std::make_optional<std::string>(db_path);
    }
    return std::nullopt;
#endif
}

MonomerDatabase::MonomerDatabase(std::string_view database_path)
{
    static constexpr std::string_view err_msg_template{
        "Problem opening database: '{}'"};

// For now, we only support non-emscripten builds.
#ifdef __EMSCRIPTEN__
    database_path = ":memory:";
#endif
    if (auto rc = sqlite3_open_v2(database_path.data(), &m_db,
                                  SQLITE_OPEN_READONLY, nullptr);
        rc != SQLITE_OK) {
        throw std::runtime_error(fmt::format(err_msg_template, database_path));
    }
}

namespace
{
template <class value_getter_t> [[nodiscard]] auto
get_info_from_database(sqlite3* db, const std::string& sql_command_template,
                       value_getter_t return_val_converter)
{
    auto sqlite3_stmt_deleter = [](auto* stmt) { sqlite3_finalize(stmt); };
    std::unique_ptr<sqlite3_stmt, decltype(sqlite3_stmt_deleter)> statement{
        nullptr, sqlite3_stmt_deleter};

    // storing this variable here instead of on the db class to improve
    // readability.
    //
    // NOTE: THE ORDER OF THE TABLE NAMES MATTERS. WE WANT TO PRIORITIZE ENTRIES
    // IN THE USER_MONOMERS TABLE OVER ENTRIES IN THE CORE_MONOMERS TABLE.
    constexpr std::array<std::string_view, 2> table_names{"user_monomers",
                                                          "core_monomers"};

    // NOTE: Return the first result you see
    for (const auto& table : table_names) {
        auto sql_command =
            fmt::format(fmt::runtime(sql_command_template), table);
        sqlite3_stmt* stmt = nullptr;
        if (auto rc =
                sqlite3_prepare_v2(db, sql_command.c_str(), -1, &stmt, NULL);
            rc != SQLITE_OK) {
            continue;
        }

        statement.reset(stmt);
        if (auto rc = sqlite3_step(stmt);
            rc != SQLITE_DONE && rc != SQLITE_OK) {
            return return_val_converter(stmt);
        };
    }

    // this is an easier way to return nullopt from this templated function
    return std::invoke_result_t<value_getter_t, sqlite3_stmt*>{};
}
} // namespace

MonomerDatabase::~MonomerDatabase()
{
    sqlite3_close_v2(m_db);
}

MonomerDatabase::string_t
MonomerDatabase::getMonomerSmiles(std::string monomer_id,
                                  ChainType monomer_type)
{
    auto get_sql_command = [&]() -> std::string {
        static constexpr std::string_view template_{
            "SELECT {0} FROM {{}} WHERE {1}='{2}' AND {3}='{4}';"};
        static constexpr std::string_view smiles_column{"SMILES"};
        static constexpr std::string_view symbol_column{"SYMBOL"};
        static constexpr std::string_view type_column{"POLYMERTYPE"};

        auto type_value = toString(monomer_type);
        boost::algorithm::trim(monomer_id);

        return fmt::format(template_, smiles_column, type_column, type_value,
                           symbol_column, monomer_id);
    };

    auto smiles_getter = [](sqlite3_stmt* stmt) -> string_t {
        auto result = sqlite3_column_text(stmt, 0);
        return std::make_optional(reinterpret_cast<const char*>(result));
    };

    return get_info_from_database(m_db, get_sql_command(), smiles_getter);
}

[[nodiscard]] MonomerDatabase::helm_info_t
MonomerDatabase::getHelmInfo(const std::string& pdb_code)
{
    auto get_sql_command = [&]() -> std::string {
        static constexpr std::string_view template_{
            "SELECT {0}, {1}, {2} FROM {{}} WHERE {3}='{4}';"};
        static constexpr std::string_view symbol_column{"SYMBOL"};
        static constexpr std::string_view type_column{"POLYMERTYPE"};
        static constexpr std::string_view smiles_column{"SMILES"};
        static constexpr std::string_view pdb_code_column{"PDBCODE"};

        return fmt::format(template_, symbol_column, type_column, smiles_column,
                           pdb_code_column, pdb_code);
    };

    auto helm_info_getter = [](sqlite3_stmt* stmt) -> helm_info_t {
        // clang-format off
        auto symbol = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        auto polymer_type = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        auto smiles = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 2));
        helm_info_t::value_type helm_info{symbol, smiles, toChainType(polymer_type)};
        return std::make_optional(helm_info);
        // clang-format on
    };

    return get_info_from_database(m_db, get_sql_command(), helm_info_getter);
}

[[nodiscard]] MonomerDatabase::string_t
MonomerDatabase::getPdbCode(const std::string& helm_symbol,
                            ChainType monomer_type)
{
    auto get_sql_command = [&]() -> std::string {
        static constexpr std::string_view template_{
            "SELECT {0} FROM {{}} WHERE {1}='{2}' AND {3}='{4}';"};
        static constexpr std::string_view pdb_code_column{"PDBCODE"};
        static constexpr std::string_view symbol_column{"SYMBOL"};
        static constexpr std::string_view type_column{"POLYMERTYPE"};

        return fmt::format(template_, pdb_code_column, symbol_column,
                           helm_symbol, type_column, toString(monomer_type));
    };

    auto pdb_code_getter = [](sqlite3_stmt* stmt) -> string_t {
        auto result = sqlite3_column_text(stmt, 0);
        return std::make_optional(reinterpret_cast<const char*>(result));
    };

    return get_info_from_database(m_db, get_sql_command(), pdb_code_getter);
}

} // namespace rdkit_extensions
} // namespace schrodinger
