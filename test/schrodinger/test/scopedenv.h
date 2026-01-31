// Set environment variable within a tight scope. Useful for test and library
// code to restore a setting of an environment variable after scope exit.
// ScopedEnvVar will retain unset or set to empty string when restoring
// enviroment.

#pragma once

#include <QByteArray>
#include <iostream>
#include <sstream>

namespace schrodinger
{
namespace test
{

class ScopedEnvVar
{
  private:
    std::string m_variable_name;
    QByteArray m_initial_value;

  public:
    // Restore setting of variable name after scope.
    ScopedEnvVar(const std::string& variable_name) :
        m_variable_name(variable_name)
    {
        m_initial_value = qgetenv(m_variable_name.c_str());
    }
    // Set variable name to new_value and restore to original value of
    // variable name on destruction of the class.
    explicit ScopedEnvVar(const std::string& variable_name,
                          const std::string& new_value) :
        m_variable_name(variable_name)
    {
        m_initial_value = qgetenv(m_variable_name.c_str());
        auto err = qputenv(m_variable_name.c_str(), new_value.c_str());
        if (err == false) {
            std::stringstream err;
            err << "Could not set env variable " << m_variable_name;
            err << " to " << new_value;
            throw std::runtime_error(err.str().c_str());
        }
    }
    ~ScopedEnvVar()
    {
        if (m_initial_value.isNull()) {
            auto err = qunsetenv(m_variable_name.c_str());
            if (err == false) {
                std::cerr << "Error calling qunsetenv for ScopedEnvVar: "
                          << m_variable_name << std::endl;
            }
        } else {
            auto err = qputenv(m_variable_name.c_str(), m_initial_value);
            if (err == false) {
                std::cerr << "Error calling qputenv for ScopedEnvVar: "
                          << m_variable_name << std::endl;
            }
        }
    }
    ScopedEnvVar(const ScopedEnvVar&) = delete;
    ScopedEnvVar& operator=(const ScopedEnvVar&) = delete;
    ScopedEnvVar(ScopedEnvVar&&) = delete;
    ScopedEnvVar& operator=(ScopedEnvVar&&) = delete;
};
// ScopedUnsetEnvVar will unset a variable_name from environment for the
// lifetime of the instance in which it is defined. It is not required for
// variable_name to be in the environment when ScopedUnsetEnvVar is defined.
class ScopedUnsetEnvVar
{
  private:
    ScopedEnvVar m_scoped_env_var;

  public:
    ScopedUnsetEnvVar(const std::string& variable_name) :
        m_scoped_env_var(variable_name)
    {
        auto err = qunsetenv(variable_name.c_str());
        if (err == false) {
            std::stringstream err;
            err << "Could not unset env variable " << variable_name;
            throw std::runtime_error(err.str().c_str());
        }
    }
    ~ScopedUnsetEnvVar(){};
    ScopedUnsetEnvVar(const ScopedUnsetEnvVar&) = delete;
    ScopedUnsetEnvVar& operator=(const ScopedUnsetEnvVar&) = delete;
    ScopedUnsetEnvVar(ScopedUnsetEnvVar&&) = delete;
    ScopedUnsetEnvVar& operator=(ScopedUnsetEnvVar&&) = delete;
};
} // namespace test
} // namespace schrodinger
