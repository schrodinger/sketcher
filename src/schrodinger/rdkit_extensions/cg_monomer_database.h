#pragma once

#include <optional>
#include <string>
#include <string_view>

#include "boost/noncopyable.hpp"

struct sqlite3;

namespace schrodinger
{
namespace rdkit_extensions
{
enum class CG_MONOMER_TYPE {
    RNA,
    PEPTIDE,

};

class [[nodiscard]] cg_monomer_database : public boost::noncopyable
{
  public:
    using monomer_smiles_t = std::optional<std::string>;

    cg_monomer_database(std::string_view database_path);

    ~cg_monomer_database();

    [[nodiscard]] monomer_smiles_t
    get_monomer_smiles(std::string monomer_id, CG_MONOMER_TYPE monomer_type);

  private:
    sqlite3* m_db;
};
} // namespace rdkit_extensions
} // namespace schrodinger
