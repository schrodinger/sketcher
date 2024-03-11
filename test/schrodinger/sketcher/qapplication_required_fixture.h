#pragma once

#include <csignal>

#include <QApplication>
#include <fmt/format.h>
#include <boost/test/unit_test.hpp>

// Use this class to construct a global fixture which would provide a
// QApplication for a whole suite of a given test file
class QApplicationRequiredFixture
{
  public:
    QApplicationRequiredFixture()
    {
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
    }

  private:
    std::unique_ptr<QApplication> d_app;
};
