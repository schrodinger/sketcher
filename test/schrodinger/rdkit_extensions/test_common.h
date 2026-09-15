#pragma once

#include <exception>
#include <filesystem>
#include <string>

#include <boost/filesystem.hpp>

#include <fmt/format.h>

#include <QByteArray> // qgetenv, qputemv, qunsetenv

#include "schrodinger/local_monomer_db_fixture.h"

BOOST_GLOBAL_FIXTURE(LocalMonomerDbFixture);

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
