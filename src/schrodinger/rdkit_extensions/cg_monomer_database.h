#pragma once

#include <optional>
#include <string>
#include <string_view>

#include "schrodinger/rdkit_extensions/definitions.h"

#include "boost/noncopyable.hpp"

struct sqlite3;

namespace schrodinger
{
namespace rdkit_extensions
{
enum class ChainType;

// Returns path of custom monomer database if it exists, otherwise returns
RDKIT_EXTENSIONS_API std::string get_custom_monomer_db_path();

std::string get_cg_monomer_db_path();

class [[nodiscard]] cg_monomer_database : public boost::noncopyable
{
  public:
    using string_t = std::optional<std::string>;
    using helm_info_t =
        std::optional<std::tuple<std::string, std::string, ChainType>>;

    cg_monomer_database(std::string_view database_path);

    ~cg_monomer_database();

    [[nodiscard]] string_t get_monomer_smiles(std::string monomer_id,
                                              ChainType monomer_type);

    [[nodiscard]] helm_info_t
    get_helm_info(const std::string& three_letter_code);

    [[nodiscard]] string_t get_pdb_code(const std::string& helm_symbol,
                                        ChainType type);

  private:
    sqlite3* m_db;
};
} // namespace rdkit_extensions
} // namespace schrodinger
