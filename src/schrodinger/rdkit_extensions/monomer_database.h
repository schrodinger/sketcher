#pragma once

#include <optional>
#include <string>
#include <string_view>
#include <vector>

#include "schrodinger/rdkit_extensions/definitions.h"

#include <boost/noncopyable.hpp>
#include <boost/filesystem/path.hpp>

struct sqlite3;

inline constexpr std::string_view CUSTOM_MONOMER_DB_PATH_ENV_VAR =
    "SCHRODINGER_CUSTOM_MONOMER_DB_PATH";

namespace schrodinger
{
namespace rdkit_extensions
{
enum class ChainType;

using opt_string_t = std::optional<std::string>;

// Returns path of custom monomer database if it exists, otherwise returns
RDKIT_EXTENSIONS_API std::optional<std::string> getMonomerDbPath();

struct RDKIT_EXTENSIONS_API MonomerInfo {
    opt_string_t symbol = std::nullopt;
    opt_string_t polymer_type = std::nullopt;
    opt_string_t natural_analog = std::nullopt;
    opt_string_t smiles = std::nullopt;
    opt_string_t core_smiles = std::nullopt;
    opt_string_t name = std::nullopt;
    opt_string_t monomer_type = std::nullopt;
    opt_string_t author = std::nullopt;
    opt_string_t pdbcode = std::nullopt;

    size_t getHash() const;
    bool areRequiredFieldsPopulated() const;
};

class RDKIT_EXTENSIONS_API MonomerDatabase : public boost::noncopyable
{
  public:
    // (HELM symbol, smiles, chain type)
    using helm_info_t =
        std::optional<std::tuple<std::string, std::string, ChainType>>;

    // Get the current DB instance. It will always contain an
    // up-to-date table with our core monomer, and MAY contain
    // an additional table with custom definitions.
    [[nodiscard]] static MonomerDatabase& instance();

    ~MonomerDatabase();

    // Populate the custom monomers table with the given SQL.
    // Note that any prior custom definitions will be dropped
    // if the SQL is successful (on failure, previous contents
    // will be preserved).
    // Also, be careful with this one, as the SQL may drop/alter
    // the table schema, which may result in a DB missing some
    // required columns!
    // @throws std::runtime_error if the SQL fails to execute,
    //    or the resulting table is missing some required columns.
    void loadMonomersFromSql(std::string_view sql);

    // Populate the custom monomers table with the definitions
    // parsed from the JSON. Note that if the JSON definitions can
    // successfully parsed, any prior custom definitions
    // will be dropped (they will be preserved if an exception is
    // thrown during the import process).
    // @throws std::runtime_error if the JSON fails to parse or
    //    if the information cannot be written into the DB.
    void loadMonomersFromJson(std::string_view json);

    // use the indicated SQLite DB file as the source for custom
    // monomer definitions.
    // @throws std::runtime_error if the DB file cannot be opened,
    //    the 'monomers' table is missing, or some required columns
    //    do not exist in the DB file.
    void loadMonomersFromSQLiteFile(boost::filesystem::path db_file);

    // Insert the custom monomers in JSON format into the custom
    // definitions table. Any previous definitions will be preserved.
    // @throws std::runtime_error if the JSON fails to parse or
    //    if the information cannot be written into the DB.
    void insertMonomersFromJson(std::string_view json);

    void resetMonomerDefinitions();

    [[nodiscard]] std::vector<std::string> getDbFields() const;

    [[nodiscard]] MonomerInfo getMonomerInfo(const std::string& monomer_id,
                                             ChainType polymer_type) const;

    [[nodiscard]] opt_string_t getMonomerSmiles(const std::string& monomer_id,
                                                ChainType polymer_type) const;

    [[nodiscard]] opt_string_t getNaturalAnalog(const std::string& monomer_id,
                                                ChainType polymer_type) const;

    [[nodiscard]] opt_string_t getPdbCode(const std::string& monomer_id,
                                          ChainType polymer_type) const;

    [[nodiscard]] helm_info_t
    getHelmInfo(const std::string& three_letter_code) const;

    // Return all information stored in the currently active DBs.
    // Note that definitions may be duplicated between the "core"
    // and the "custom" DBs.
    // @param include_core : whether to include the core monomers
    //   in the output or just return the custom ones. Default: true.
    // @returns a JSON string with the content of the DBs.
    [[nodiscard]] std::string
    getMonomerDefinitions(bool include_core = true) const;

    // Regenerate the SMILES and CORE_SMILES columns in the DBs.
    // This should be used after loading SQL/JSON/a SQlite file
    // to make sure the SMILES strings are correctly canonicalized.
    void canonicalizeSmilesFields(bool include_core = false);

    [[nodiscard]] std::vector<std::pair<std::string, std::string>>
    getAllSMILES() const;

  private:
    // this is private because we don't want to allow managing
    // just any db -- we require the proper schema!
    MonomerDatabase();

    // swaps the SQLite DB we use to read the custom monomers
    // from. The instance will take ownership of the DB and close
    // the current DB instance if there's one.
    void swap_custom_monomers_db(sqlite3* db);

    // Checks whether the db has the schema that is required
    // by this class (i.e. the table(s) and columns a the
    // core monomers db)
    [[nodiscard]] std::vector<std::string> check_db(sqlite3* db) const;

    void dumpToFile(sqlite3* db, boost::filesystem::path db_file) const;

    sqlite3* m_core_monomers_db = nullptr;
    sqlite3* m_custom_monomers_db = nullptr;
};
} // namespace rdkit_extensions
} // namespace schrodinger
