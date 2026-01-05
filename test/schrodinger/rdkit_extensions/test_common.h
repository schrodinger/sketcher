#pragma once

#include <exception>
#include <filesystem>
#include <string>

#include <boost/filesystem.hpp>

#include <fmt/format.h>

#include <QByteArray> // qgetenv, qputemv, qunsetenv

#include "schrodinger/rdkit_extensions/monomer_database.h"

/**
 * Gets full path to a file in the testfiles directory in the source
 * @param filename file found in the testfiles folder
 * @return full filesystem path to that file
 */
std::string testfile_path(const std::string& filename)
{
    auto path = std::filesystem::path(std::getenv("SKETCHER_SOURCE_DIR")) /
                "test" / "testfiles" / filename;
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("File not found: " +
                                 std::filesystem::absolute(path).string());
    }
    return path.string();
}

// Use a custom monomer DB file in the test directory to avoid
// creating one under ~/.schrodinger

class LocalMonomerDbFixture
{
  public:
    LocalMonomerDbFixture()
    {
        using schrodinger::rdkit_extensions::getMonomerDbPath;

        // We want to override the default path for the custom monomer db,
        // so that we don't accidentally use anything the user has defined
        // in whatever tests we run, while we still want to be able to
        // have custom definitions in the tests.
        m_custom_db = qgetenv(CUSTOM_MONOMER_DB_PATH_ENV_VAR.data());

        auto test_dir = boost::filesystem::current_path();
        auto db_fname =
            boost::unit_test::framework::master_test_suite().p_name.value +
            "_custom_monomers.db";

        // Make sure the  file doesn't exist, so that we aren't using
        // preexisting data in the tests
        auto tmp_db = test_dir / db_fname;
        if (boost::filesystem::exists(tmp_db)) {
            boost::filesystem::remove(tmp_db);
        }

        // We want to exit early if we can't set up
        // CUSTOM_MONOMER_DB_PATH_ENV_VAR because we don't want to accidentally
        // overwrite the user's custom monomer DB.
        // For whatever reasons, errors inside a global fixture are not
        // reported as a test failure (though they abort the test run
        // and exit with status code 200).

        // We want to keep this one while the fixture is in effect,
        // so that we don't pull the rug from under the env var.
        test_custom_monomer_db = tmp_db.string();

        // We ust Qt for env vars because putenv/setenv are not portable
        // (they don't exist in MSVC)
        if (qputenv(CUSTOM_MONOMER_DB_PATH_ENV_VAR.data(),
                    test_custom_monomer_db.c_str()) == 0) {
            throw std::runtime_error("\n\nError: Failed setting temporary "
                                     "custom monomer db file.\n\n");
        }

        if (auto check = getMonomerDbPath(); !check.has_value()) {
            throw std::runtime_error(
                "\n\nError: getMonomerDbPath() did not return a value.\n\n");
        } else if (*check != test_custom_monomer_db) {
            auto msg = fmt::format("\n\nError: getMonomerDbPath() does not "
                                   "match the env var: {} != {}.\n\n",
                                   *check, test_custom_monomer_db);
            throw std::runtime_error(msg);
        }
    }

    ~LocalMonomerDbFixture()
    {
        // Restore whatever default we had before the test,
        // or at least clean up what this fixture set up
        if (m_custom_db.isNull()) {
            qunsetenv(CUSTOM_MONOMER_DB_PATH_ENV_VAR.data());
        } else {
            qputenv(CUSTOM_MONOMER_DB_PATH_ENV_VAR.data(), m_custom_db);
        }
    }

    static std::string test_custom_monomer_db;

  private:
    QByteArray m_custom_db;
};

std::string LocalMonomerDbFixture::test_custom_monomer_db;

BOOST_GLOBAL_FIXTURE(LocalMonomerDbFixture);