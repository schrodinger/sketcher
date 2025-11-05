/* -------------------------------------------------------------------------
 * Common schrodinger::test:: fixtures for boost unit tests
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#include <iostream>
#include <sstream>
#include <streambuf>

namespace schrodinger
{
namespace test
{

/**
 * @brief Temporarily silence standard output and error streams.
 *
 * The `silence_stdlog` fixture captures all writes to `std::cout` and
 * `std::cerr` during its lifetime by redirecting them to internal
 * stringstreams. When the fixture is destroyed, the original output and error
 * stream buffers are restored.
 *
 * Example usage:
 * @code
 * // For an individual test case:
 * BOOST_FIXTURE_TEST_CASE(test_name, silence_stdlog) {
 *     // ... test logic ...
 * }
 *
 * // To silence all tests in a module:
 * BOOST_GLOBAL_FIXTURE(silence_stdlog)
 * @endcode
 */
struct silence_stdlog {
    /**
     * @brief Constructs the fixture and redirects stdout and stderr.
     *
     * Saves the current `std::cout` and `std::cerr` stream buffers, then
     * replaces them with internal stringstreams to suppress or capture output.
     */
    silence_stdlog() :
        m_prev_errbuf(std::cerr.rdbuf()),
        m_prev_outbuf(std::cout.rdbuf())
    {
        std::cerr.rdbuf(m_errstream.rdbuf());
        std::cout.rdbuf(m_outstream.rdbuf());
    }

    /**
     * @brief Restores the original stdout and stderr stream buffers.
     *
     * Ensures that any redirection done in the constructor is reverted.
     * Clears the error state of `std::cout` and `std::cerr` after restoration.
     */
    ~silence_stdlog()
    {
        std::cerr.rdbuf(m_prev_errbuf);
        std::cout.rdbuf(m_prev_outbuf);

        std::cerr.clear();
        std::cout.clear();
    }

    /** @brief Pointer to the original `std::cerr` stream buffer. */
    std::streambuf* m_prev_errbuf;

    /** @brief Pointer to the original `std::cout` stream buffer. */
    std::streambuf* m_prev_outbuf;

    /** @brief Internal stream used to capture redirected stderr output. */
    std::stringstream m_errstream;

    /** @brief Internal stream used to capture redirected stdout output. */
    std::stringstream m_outstream;
};

} // namespace test
} // namespace schrodinger
