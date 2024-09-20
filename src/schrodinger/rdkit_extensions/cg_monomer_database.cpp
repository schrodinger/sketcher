#include "schrodinger/rdkit_extensions/cg_monomer_database.h"

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
            "SELECT {0} FROM {1} WHERE {2}='{3}' AND {4}='{5}';"};
        static constexpr std::string_view core_monomers_table{"core_monomers"};
        static constexpr std::string_view smiles_column{"SMILES"};
        static constexpr std::string_view symbol_column{"SYMBOL"};
        static constexpr std::string_view type_column{"POLYMERTYPE"};

        auto type_value = to_string(monomer_type);
        boost::algorithm::trim(monomer_id);

        return fmt::format(template_, smiles_column, core_monomers_table,
                           type_column, type_value, symbol_column, monomer_id);
    };

    auto sqlite3_stmt_deleter = [](auto* stmt) { sqlite3_finalize(stmt); };
    std::unique_ptr<sqlite3_stmt, decltype(sqlite3_stmt_deleter)> statement{
        nullptr, sqlite3_stmt_deleter};

    sqlite3_stmt* stmt = nullptr;
    if (auto rc = sqlite3_prepare_v2(m_db, get_sql_command().c_str(), -1, &stmt,
                                     NULL);
        rc != SQLITE_OK) {
        return std::nullopt;
    }

    statement.reset(stmt);
    if (auto rc = sqlite3_step(stmt); rc != SQLITE_DONE && rc != SQLITE_OK) {
        auto get_column_text = [&stmt](int idx) -> std::string {
            auto result = sqlite3_column_text(stmt, idx);
            return {reinterpret_cast<const char*>(result)};
        };

        return std::make_optional(get_column_text(0));
    }
    return std::nullopt;
}

[[nodiscard]] cg_monomer_database::helm_info_t
cg_monomer_database::get_helm_info(const std::string& pdb_code)
{
    auto get_sql_command = [&]() -> std::string {
        static constexpr std::string_view template_{
            "SELECT {0}, {1} FROM {2} WHERE {3}='{4}';"};
        static constexpr std::string_view core_monomers_table{"core_monomers"};
        static constexpr std::string_view symbol_column{"SYMBOL"};
        static constexpr std::string_view type_column{"POLYMERTYPE"};
        static constexpr std::string_view pdb_code_column{"PDBCODE"};

        return fmt::format(template_, symbol_column, type_column,
                           core_monomers_table, pdb_code_column, pdb_code);
    };

    auto sqlite3_stmt_deleter = [](auto* stmt) { sqlite3_finalize(stmt); };
    std::unique_ptr<sqlite3_stmt, decltype(sqlite3_stmt_deleter)> statement{
        nullptr, sqlite3_stmt_deleter};

    sqlite3_stmt* stmt = nullptr;
    if (auto rc = sqlite3_prepare_v2(m_db, get_sql_command().c_str(), -1, &stmt,
                                     NULL);
        rc != SQLITE_OK) {
        return std::nullopt;
    }

    statement.reset(stmt);
    if (auto rc = sqlite3_step(stmt); rc != SQLITE_DONE && rc != SQLITE_OK) {
        auto get_column_text = [&stmt](int idx) -> std::string {
            auto result = sqlite3_column_text(stmt, idx);
            return {reinterpret_cast<const char*>(result)};
        };
        auto monomer_type = to_chain_type(get_column_text(1));
        std::pair<std::string, ChainType> helm_info = {
            std::string(get_column_text(0)), monomer_type};
        return std::make_optional(helm_info);
    }
    return std::nullopt;
}

} // namespace rdkit_extensions
} // namespace schrodinger
