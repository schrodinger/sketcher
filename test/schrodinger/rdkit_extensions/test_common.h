#pragma once

#include <cstdlib>
#include <exception>
#include <string>

#include <filesystem>

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

/**
 * Set the default monomer db path
 */
void set_default_monomer_db_path()
{
    auto path = std::filesystem::path(std::getenv("SKETCHER_SOURCE_DIR")) /
                "data" / "core_monomerlib.db";
    if (!std::filesystem::exists(path)) {
        throw std::runtime_error("Monomer database not found: " +
                                 std::filesystem::absolute(path).string());
    }

    // Set the environment variable
    auto name = "SCHRODINGER_DEFAULT_MONOMER_DB_PATH";
    auto value = path.string();
#ifdef _WIN32
    _putenv_s(name, value.c_str());
#else
    int overwrite = 1;
    setenv(name, value.c_str(), overwrite);
#endif
}