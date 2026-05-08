/* -------------------------------------------------------------------------
 * Utility to convert SMILES to HELM format
 *
 * Reads SMILES strings from stdin (one per line, optionally followed by a
 * space and a name), converts them to monomeric representation, and outputs
 * HELM format to stdout (with optional name).
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/filesystem.hpp>

#include <QByteArray>
#include <QtEnvironmentVariables> // qputenv

#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/RWMol.h>

#include "schrodinger/rdkit_extensions/atomistic_conversions.h"
#include "schrodinger/rdkit_extensions/convert.h"
#include "schrodinger/rdkit_extensions/file_format.h"
#include "schrodinger/rdkit_extensions/monomer_database.h"

using namespace schrodinger::rdkit_extensions;

namespace
{
// Global variable to track temporary database file for cleanup
boost::filesystem::path g_temp_db_path;

void cleanup_temp_db()
{
    if (!g_temp_db_path.empty() && boost::filesystem::exists(g_temp_db_path)) {
        boost::filesystem::remove(g_temp_db_path);
    }
}

struct Options {
    std::string db_path;
};

auto help_message = R"(
Convert SMILES to HELM format.

Reads SMILES strings from stdin (one per line, optionally followed by a space
and a name), converts them to monomeric representation, and outputs HELM
format to stdout (with optional name).

Options:
  --db FILE     Load custom monomer database from FILE (.sqlite/.db or .json)
  -h, --help    Show this help message and exit\n";
)";

void print_usage(const char* program_name)
{
    std::cout << "Usage: " << program_name << " [OPTIONS]\n" << help_message;
}

// Parse command-line arguments
Options parse_args(int argc, char* argv[])
{
    Options options;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];

        if (arg == "-h" || arg == "--help") {
            print_usage(argv[0]);
            std::exit(0);
        } else if (arg == "--db") {
            if (i + 1 >= argc) {
                std::cerr << "Error: --db requires a file path\n";
                print_usage(argv[0]);
                std::exit(1);
            }
            options.db_path = argv[++i];
        } else {
            std::cerr << "Error: Unknown option: " << arg << "\n";
            print_usage(argv[0]);
            std::exit(1);
        }
    }

    return options;
}

// Load a JSON or SQLite monomer database into the MonomerDatabase singleton.
void load_custom_database(const std::string& db_path)
{
    namespace fs = boost::filesystem;

    fs::path path(db_path);
    if (!fs::exists(path)) {
        throw std::runtime_error("Database file does not exist: " + db_path);
    }

    auto& db = MonomerDatabase::instance();
    std::string ext = path.extension().string();

    if (ext == ".sqlite" || ext == ".db") {
        db.loadMonomersFromSQLiteFile(path);
    } else if (ext == ".json") {
        // For JSON files, set environment variable to a temporary .db file
        // before loading, so the database gets written there instead of the
        // default location
        g_temp_db_path =
            fs::temp_directory_path() /
            fs::unique_path("smiles_to_helm_%%%%-%%%%-%%%%-%%%%.db");

        // Register cleanup handler
        std::atexit(cleanup_temp_db);

        // Use Qt for env vars because putenv/setenv are not portable
        qputenv(CUSTOM_MONOMER_DB_PATH_ENV_VAR.data(),
                g_temp_db_path.string().c_str());

        // Read JSON file
        std::ifstream file(db_path);
        if (!file.is_open()) {
            throw std::runtime_error("Failed to open file: " + db_path);
        }
        std::stringstream buffer;
        buffer << file.rdbuf();

        // Load JSON into database (this will write to g_temp_db_path)
        db.loadMonomersFromJson(buffer.str());
    } else {
        throw std::runtime_error("Unsupported file extension: " + ext +
                                 " (expected .sqlite or .json)");
    }
}
} // namespace

int main(int argc, char* argv[])
{
    auto options = parse_args(argc, argv);

    try {
        // Load custom database if provided
        if (!options.db_path.empty()) {
            load_custom_database(options.db_path);
        }

    } catch (const std::exception& e) {
        std::cerr << "Error loading database: " << e.what() << std::endl;
        return 1;
    }

    std::string line;
    auto status = 0;

    while (std::getline(std::cin, line)) {
        if (line.empty()) {
            continue;
        }

        try {
            auto mol = to_rdkit(line, Format::SMILES);
            auto name = std::string{};
            mol->getPropIfPresent("_Name", name);

            const RDKit::ROMol& mol_ref = *mol;
            auto monomer_mol = toMonomeric(mol_ref);

            std::string helm = to_string(*monomer_mol, Format::HELM);
            std::cout << helm;
            if (!name.empty()) {
                std::cout << " " << name;
            }
            std::cout << std::endl;

        } catch (const std::exception& e) {
            std::cerr << "Error processing SMILES '" << line
                      << "': " << e.what() << std::endl;
            // Keep going, but exit with error status when done.
            status = 1;
        }
    }

    return status;
}
