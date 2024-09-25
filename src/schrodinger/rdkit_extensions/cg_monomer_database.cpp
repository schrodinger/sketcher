#include "schrodinger/rdkit_extensions/cg_monomer_database.h"

#include <array>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include "boost/algorithm/string/trim.hpp"
#include "boost/noncopyable.hpp"

#include "schrodinger/rdkit_extensions/coarse_grain.h" // ChainType

#ifndef __EMSCRIPTEN__
#include "mmfile.h"
#include "schrodinger/path.h"
#endif

#include <fmt/format.h>

#include "sqlite3.h"

using namespace std::string_view_literals;

namespace schrodinger
{
namespace rdkit_extensions
{

std::string get_cg_monomer_db_path()
{
#ifdef __EMSCRIPTEN__
    throw std::logic_error(
        "CG Monomers are not yet supported in WASM Sketcher");
#else
    // First check for custom monomer database in .schrodinger, then fall back
    // to default
    auto custom_db_path = boost::filesystem::path(mmfile_get_directory_path(
                              DirectoryName::MMFILE_LOCAL_APPDATA)) /
                          "helm/custom_monomerlib.db";
    if (boost::filesystem::exists(custom_db_path)) {
        return custom_db_path.string();
    }

    auto db_path =
        path::product_dir("mmshare") / "data/helm/core_monomerlib.db";
    return db_path.string();
#endif
}

cg_monomer_database::cg_monomer_database(std::string_view database_path)
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

cg_monomer_database::~cg_monomer_database()
{
    sqlite3_close_v2(m_db);
}

cg_monomer_database::monomer_smiles_t
cg_monomer_database::get_monomer_smiles(std::string monomer_id,
                                        ChainType monomer_type)
{
    auto get_sql_command = [&]() -> std::string {
        static constexpr std::string_view template_{
            "SELECT {0} FROM {{}} WHERE {1}='{2}' AND {3}='{4}';"};
        static constexpr std::string_view smiles_column{"SMILES"};
        static constexpr std::string_view symbol_column{"SYMBOL"};
        static constexpr std::string_view type_column{"POLYMERTYPE"};

        auto type_value = to_string(monomer_type);
        boost::algorithm::trim(monomer_id);

        return fmt::format(template_, smiles_column, type_column, type_value,
                           symbol_column, monomer_id);
    };

    auto smiles_getter = [](sqlite3_stmt* stmt) -> monomer_smiles_t {
        auto result = sqlite3_column_text(stmt, 0);
        return std::make_optional(reinterpret_cast<const char*>(result));
    };

    return get_info_from_database(m_db, get_sql_command(), smiles_getter);
}

[[nodiscard]] cg_monomer_database::helm_info_t
cg_monomer_database::get_helm_info(const std::string& pdb_code)
{
    auto get_sql_command = [&]() -> std::string {
        static constexpr std::string_view template_{
            "SELECT {0}, {1} FROM {{}} WHERE {2}='{3}';"};
        static constexpr std::string_view symbol_column{"SYMBOL"};
        static constexpr std::string_view type_column{"POLYMERTYPE"};
        static constexpr std::string_view pdb_code_column{"PDBCODE"};

        return fmt::format(template_, symbol_column, type_column,
                           pdb_code_column, pdb_code);
    };

    auto helm_info_getter = [](sqlite3_stmt* stmt) -> helm_info_t {
        // clang-format off
        auto symbol = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 0));
        auto polymer_type = reinterpret_cast<const char*>(sqlite3_column_text(stmt, 1));
        helm_info_t::value_type helm_info{symbol, to_chain_type(polymer_type)};
        return std::make_optional(helm_info);
        // clang-format on
    };

    return get_info_from_database(m_db, get_sql_command(), helm_info_getter);
}

} // namespace rdkit_extensions
} // namespace schrodinger
