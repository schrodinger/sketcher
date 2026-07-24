// @copyright Schrodinger, LLC - All Rights Reserved

#pragma once

#include <csignal>
#include <vector>

#include <QApplication>
#include <QByteArray>
#include <fmt/format.h>
#include <boost/test/unit_test.hpp>
#include <cstdlib>

#include "schrodinger/sketcher/font_loader.h"

/// @return true if there is a display
static bool has_display()
{
#ifdef __linux__
    return getenv("DISPLAY") != nullptr;
#elif __APPLE__
    return getenv("SSH_CONNECTION") == nullptr;
#else
    // Always DISPLAY available on Windows.
    return true;
#endif
}

// Use this class to construct a global fixture which would provide a
// QApplication for a whole suite of a given test file
class QApplicationRequiredFixture
{
  public:
    QApplicationRequiredFixture()
    {
        if (!has_display()) {
            BOOST_TEST_MESSAGE("Skipping tests that require a display.");
            exit(0);
        }

        if (QCoreApplication::instance() != nullptr) {
            throw std::runtime_error(
                "QApplication instance already exists, multiple QApplication "
                "instances are not allowed.");
        }

        // assert() raises a SIGABRT signal, which is not an exception
        // and bypasses the destructor, so add a handler for it.
        signal(SIGABRT, [](int signal) {
            auto test_name =
                boost::unit_test::framework::current_test_case().p_name.get();
            throw std::runtime_error(fmt::format(
                "A SIGABRT signal was raised from '{}'", test_name));
        });

        // Qt requires argc and argv to stay valid for the entire lifetime of
        // the QApplication object.
        const auto& test_suite =
            boost::unit_test::framework::master_test_suite();
        d_arguments.reserve(test_suite.argc + 2);
        for (int i = 0; i < test_suite.argc; ++i) {
            d_arguments.emplace_back(test_suite.argv[i]);
        }
        d_arguments.emplace_back("--platform");
        d_arguments.emplace_back("offscreen");

        d_argv.reserve(d_arguments.size());
        for (auto& argument : d_arguments) {
            d_argv.push_back(argument.data());
        }
        d_argc = static_cast<int>(d_argv.size());
        d_app.reset(new QApplication(d_argc, d_argv.data()));

        // Initialize Qt resources (fonts, icons, etc.)
#ifdef SKETCHER_STATIC_DEFINE
        Q_INIT_RESOURCE(sketcher);
#endif
        // load the Arimo font so that font width calculations result in
        // expected values
        schrodinger::sketcher::load_font_resources();
    }

  private:
    int d_argc = 0;
    std::vector<QByteArray> d_arguments;
    std::vector<char*> d_argv;
    std::unique_ptr<QApplication> d_app;
};
