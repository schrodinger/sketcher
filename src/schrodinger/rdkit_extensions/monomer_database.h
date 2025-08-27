#pragma once

#include <optional>
#include <string>
#include <string_view>

#include "schrodinger/rdkit_extensions/definitions.h"

#include "boost/noncopyable.hpp"

struct sqlite3;

const std::string CUSTOM_MONOMER_DB_PATH_ENV_VAR =
    "SCHRODINGER_CUSTOM_MONOMER_DB_PATH";
const std::string DEFAULT_MONOMER_DB_PATH_ENV_VAR =
    "SCHRODINGER_DEFAULT_MONOMER_DB_PATH";

namespace schrodinger
{
namespace rdkit_extensions
{
enum class ChainType;

// Returns path of custom monomer database if it exists, otherwise returns
RDKIT_EXTENSIONS_API std::optional<std::string> getCustomMonomerDbPath();

RDKIT_EXTENSIONS_API std::optional<std::string> getMonomerDbPath();

class RDKIT_EXTENSIONS_API MonomerDatabase : public boost::noncopyable
{
  public:
    using string_t = std::optional<std::string>;
    // (HELM symbol, smiles, chain type)
    using helm_info_t =
        std::optional<std::tuple<std::string, std::string, ChainType>>;

    MonomerDatabase(std::string_view database_path);

    ~MonomerDatabase();

    [[nodiscard]] string_t getMonomerSmiles(std::string monomer_id,
                                            ChainType monomer_type) const;

    [[nodiscard]] std::string getNaturalAnalog(std::string monomer_id,
                                               ChainType monomer_type) const;

    [[nodiscard]] helm_info_t
    getHelmInfo(const std::string& three_letter_code) const;

    [[nodiscard]] string_t getPdbCode(const std::string& helm_symbol,
                                      ChainType type) const;

  private:
    sqlite3* m_db;
};
} // namespace rdkit_extensions
} // namespace schrodinger
