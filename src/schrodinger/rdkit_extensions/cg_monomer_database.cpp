#include "schrodinger/rdkit_extensions/cg_monomer_database.h"

#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>

#include "boost/algorithm/string/trim.hpp"
#include "boost/noncopyable.hpp"

#include <fmt/format.h>

#include "sqlite3.h"

using namespace std::string_view_literals;

namespace schrodinger
{
namespace rdkit_extensions
{

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
                                        CG_MONOMER_TYPE monomer_type)
{
    auto get_sql_command = [&]() -> std::string {
        static constexpr std::string_view template_{
            "SELECT {0} FROM {1} WHERE {2}='{3}' AND {4}='{5}';"};
        static constexpr std::string_view core_monomers_table{"core_monomers"};
        static constexpr std::string_view smiles_column{"SMILES"};
        static constexpr std::string_view symbol_column{"SYMBOL"};
        static constexpr std::string_view type_column{"POLYMERTYPE"};

        auto type_value =
            (monomer_type == CG_MONOMER_TYPE::RNA ? "RNA"sv : "PEPTIDE"sv);
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

} // namespace rdkit_extensions
} // namespace schrodinger
