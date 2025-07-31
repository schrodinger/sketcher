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
    std::filesystem::path path;
    if (auto src_path = std::getenv("SKETCHER_SOURCE_DIR")) {
        path =
            std::filesystem::path(src_path) / "test" / "testfiles" / filename;
    }
    if (auto src_path = std::getenv("SCHRODINGER_SRC")) {
        path = std::filesystem::path(src_path) / "mmshare" / "test" /
               "testfiles" / "rdkit_extensions" / filename;
    }
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
    std::filesystem::path path;
    if (auto src_path = std::getenv("SKETCHER_SOURCE_DIR")) {
        path = std::filesystem::path(src_path) / "data" / "core_monomerlib.db";
    }
    if (auto src_path = std::getenv("SCHRODINGER_SRC")) {
        path = std::filesystem::path(src_path) / "mmshare" / "data" / "helm" /
               "core_monomerlib.db";
    }
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