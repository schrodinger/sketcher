#pragma once

#include <stdexcept>
#include <string>

#include <boost/test/unit_test.hpp>
#include <fmt/format.h>
#include <QByteArray>
#include <QTemporaryDir>

#include "schrodinger/rdkit_extensions/monomer_database.h"
// Give each test process a unique database, even concurrent runs of the same
// executable. Never read or overwrite the user's custom monomer database.

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

        if (!m_temp_dir.isValid()) {
            throw std::runtime_error(
                "Failed creating temporary monomer DB directory.");
        }
        // We want to exit early if we can't set up
        // CUSTOM_MONOMER_DB_PATH_ENV_VAR because we don't want to accidentally
        // overwrite the user's custom monomer DB.
        // For whatever reasons, errors inside a global fixture are not
        // reported as a test failure (though they abort the test run
        // and exit with status code 200).

        // We want to keep this one while the fixture is in effect,
        // so that we don't pull the rug from under the env var.
        test_custom_monomer_db =
            m_temp_dir.filePath("custom_monomers.db").toStdString();

        // We use Qt for env vars because putenv/setenv are not portable
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
        // Close file-backed connections before removing the temporary
        // directory.
        schrodinger::rdkit_extensions::MonomerDatabase::instance()
            .resetMonomerDefinitions();

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
    QTemporaryDir m_temp_dir;
    QByteArray m_custom_db;
};

std::string LocalMonomerDbFixture::test_custom_monomer_db;
