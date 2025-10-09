// @copyright Schrodinger, LLC - All Rights Reserved

#pragma once

#include <csignal>

#include <QApplication>
#include <fmt/format.h>
#include <boost/test/unit_test.hpp>
#include <cstdlib>

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
        // the QApplication object, so pass boost master test suite variables.
        d_app.reset(new QApplication(
            boost::unit_test::framework::master_test_suite().argc,
            boost::unit_test::framework::master_test_suite().argv));

        // Initialize Qt resources (fonts, icons, etc.)
#ifdef SKETCHER_STATIC_DEFINE
        Q_INIT_RESOURCE(sketcher);
#endif
    }

  private:
    std::unique_ptr<QApplication> d_app;
};
