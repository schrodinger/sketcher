/* -------------------------------------------------------------------------
 * Common schrodinger::test:: utility to check exception messages
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <string>
#include <sstream>
#include <regex>
#include <boost/test/unit_test.hpp>

/**
 * Public definitions
 * see /mmshare/test/schrodinger/test/test_checkexceptionmsg.cpp for an example
 * on how to use
 */

#define TEST_WARN_EXCEPTION_MSG_SUBSTR(S, E, M)  \
    BOOST_WARN_EXCEPTION(                        \
        S, E,                                    \
        schrodinger::test::CheckExceptionMsg<E>( \
            M, schrodinger::test::CheckExceptionMsg<E>::TestSeverity::WARN));

#define TEST_CHECK_EXCEPTION_MSG_SUBSTR(S, E, M) \
    BOOST_CHECK_EXCEPTION(                       \
        S, E,                                    \
        schrodinger::test::CheckExceptionMsg<E>( \
            M, schrodinger::test::CheckExceptionMsg<E>::TestSeverity::CHECK));

#define TEST_REQUIRE_EXCEPTION_MSG_SUBSTR(S, E, M) \
    BOOST_REQUIRE_EXCEPTION(                       \
        S, E,                                      \
        schrodinger::test::CheckExceptionMsg<E>(   \
            M,                                     \
            schrodinger::test::CheckExceptionMsg<E>::TestSeverity::REQUIRE));

#define TEST_WARN_EXCEPTION_MSG_REGEX(S, E, M)                             \
    BOOST_WARN_EXCEPTION(S, E,                                             \
                         schrodinger::test::CheckExceptionMsgRegex<E>(     \
                             M, schrodinger::test::CheckExceptionMsgRegex< \
                                    E>::TestSeverity::WARN));

#define TEST_CHECK_EXCEPTION_MSG_REGEX(S, E, M)                             \
    BOOST_CHECK_EXCEPTION(S, E,                                             \
                          schrodinger::test::CheckExceptionMsgRegex<E>(     \
                              M, schrodinger::test::CheckExceptionMsgRegex< \
                                     E>::TestSeverity::CHECK));

#define TEST_REQUIRE_EXCEPTION_MSG_REGEX(S, E, M)                             \
    BOOST_REQUIRE_EXCEPTION(S, E,                                             \
                            schrodinger::test::CheckExceptionMsgRegex<E>(     \
                                M, schrodinger::test::CheckExceptionMsgRegex< \
                                       E>::TestSeverity::REQUIRE));

namespace schrodinger
{
namespace test
{

template <class T> class CheckExceptionMsg
{
  public:
    enum class TestSeverity { WARN, CHECK, REQUIRE };

    CheckExceptionMsg(std::string msg, TestSeverity severity) :
        m_severity(severity),
        m_expected(std::move(msg)),
        m_match(false)
    {
    }

    ~CheckExceptionMsg()
    {
        // This message should be printer after BOOST_REQUIRE_EXCEPTION's
        std::stringstream ss;
        ss << '\"' << m_expected << "\" not found in \"" << m_message << '\"';

        switch (m_severity) {
            case TestSeverity::WARN:
                BOOST_WARN_MESSAGE(m_match, ss.str());
                break;
            case TestSeverity::CHECK:
                BOOST_CHECK_MESSAGE(m_match, ss.str());
                break;
            case TestSeverity::REQUIRE:
                BOOST_REQUIRE_MESSAGE(m_match, ss.str());
                break;
        }
    }

    bool operator()(const T& exc)
    {
        m_message = exc.what();
        m_match = m_message.find(m_expected) != std::string::npos;

        return m_match;
    }

  private:
    const TestSeverity m_severity;
    const std::string m_expected;
    std::string m_message;
    bool m_match;
};

template <class T> class CheckExceptionMsgRegex
{
  public:
    enum class TestSeverity { WARN, CHECK, REQUIRE };

    CheckExceptionMsgRegex(std::string msg, TestSeverity severity) :
        m_severity(severity),
        m_expected(std::move(msg)),
        m_match(false)
    {
    }

    ~CheckExceptionMsgRegex()
    {
        // This message should be printer after BOOST_REQUIRE_EXCEPTION's
        std::stringstream ss;
        ss << '\"' << m_expected << "\" doesn't match \"" << m_message << '\"';

        switch (m_severity) {
            case TestSeverity::WARN:
                BOOST_WARN_MESSAGE(m_match, ss.str());
                break;
            case TestSeverity::CHECK:
                BOOST_CHECK_MESSAGE(m_match, ss.str());
                break;
            case TestSeverity::REQUIRE:
                BOOST_REQUIRE_MESSAGE(m_match, ss.str());
                break;
        }
    }

    bool operator()(const T& exc)
    {
        m_message = exc.what();
        std::smatch match;
        // usually you want to compile the re for speed. But here we actually
        // do want the re string because reporting legibly is more important
        // than performance
        std::regex expr{m_expected};
        m_match = std::regex_search(m_message, match, expr);

        return m_match;
    }

  private:
    const TestSeverity m_severity;
    const std::string m_expected;
    std::string m_message;
    bool m_match;
};

} // namespace test
} // namespace schrodinger