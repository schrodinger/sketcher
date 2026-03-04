/* -------------------------------------------------------------------------
 * Crash handler for the standalone schrodinger_sketcher application.
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#ifndef __EMSCRIPTEN__

#include "crash_handler.h"

#include <chrono>
#include <csignal>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>

#include <boost/stacktrace.hpp>
#include <fmt/format.h>

#include "schrodinger/sketcher/version.h"

#ifdef _WIN32
#include <windows.h>
#endif

namespace schrodinger
{
namespace
{

/// Get the current local time as a std::tm.
std::tm get_local_time()
{
    auto now =
        std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    return *std::localtime(&now);
}

/**
 * Return the path where this crash report should be written. Uses the system
 * temp directory with a timestamped filename.
 */
std::filesystem::path get_crash_report_path(const std::tm& tm)
{
    auto filename =
        fmt::format("sketcher_crash_{:04d}{:02d}{:02d}_{:02d}{:02d}{:02d}.log",
                    tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour,
                    tm.tm_min, tm.tm_sec);

    return std::filesystem::temp_directory_path() / filename;
}

/// Write the crash report to stderr and to a file.
void write_crash_report(std::string_view reason,
                        const boost::stacktrace::stacktrace& trace)
{
    auto tm = get_local_time();
    auto trace_str = boost::stacktrace::to_string(trace);
    auto report = fmt::format(
        "=== Schrodinger 2D Sketcher Crash Report ===\n"
        "Version: {} (build {})\n"
        "Time: {:04d}-{:02d}-{:02d} {:02d}:{:02d}:{:02d}\n"
        "Reason: {}\n"
        "\n"
        "Stack trace:\n"
        "{}\n",
        std::string(sketcher::SKETCHER_RELEASE), sketcher::SKETCHER_BUILD,
        tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min,
        tm.tm_sec, reason, trace_str);

    // Write to stderr (best-effort)
    std::cerr << report << std::flush;

    // Write to file (best-effort)
    try {
        auto path = get_crash_report_path(tm);
        std::ofstream file(path);
        if (file.is_open()) {
            file << report << std::flush;
            std::cerr << "\nCrash report written to: " << path << "\n";
        }
    } catch (...) {
        // If file writing fails, we already wrote to stderr
    }
}

/// Signal handler for fatal signals (SIGSEGV, SIGABRT, etc.)
void signal_handler(int signal_number)
{
    // Reset to default handler to avoid infinite recursion if the
    // crash report code itself triggers the same signal.
    std::signal(signal_number, SIG_DFL);

    const char* signal_name = "Unknown signal";
    switch (signal_number) {
        case SIGSEGV:
            signal_name = "SIGSEGV (Segmentation fault)";
            break;
        case SIGABRT:
            signal_name = "SIGABRT (Abort)";
            break;
        case SIGFPE:
            signal_name = "SIGFPE (Floating point exception)";
            break;
        case SIGILL:
            signal_name = "SIGILL (Illegal instruction)";
            break;
#ifndef _WIN32
        case SIGBUS:
            signal_name = "SIGBUS (Bus error)";
            break;
#endif
    }

    auto reason =
        fmt::format("Fatal signal: {} ({})", signal_name, signal_number);
    write_crash_report(reason, boost::stacktrace::stacktrace());

    // Re-raise so the OS default handler runs (core dump, WER dialog, etc.)
    std::raise(signal_number);
}

/// Terminate handler for unhandled C++ exceptions.
void terminate_handler()
{
    std::string reason = "std::terminate() called";

    // Try to extract information about the current exception
    if (auto eptr = std::current_exception()) {
        try {
            std::rethrow_exception(eptr);
        } catch (const std::exception& e) {
            reason = fmt::format("Unhandled exception: {}", e.what());
        } catch (...) {
            reason = "Unhandled exception (unknown type)";
        }
    }

    write_crash_report(reason, boost::stacktrace::stacktrace());
    std::abort();
}

#ifdef _WIN32
/**
 * Windows Structured Exception Handling (SEH) for access violations etc. that
 * bypass C++ exception handling.
 */
LONG WINAPI unhandled_exception_filter(EXCEPTION_POINTERS* exception_info)
{
    DWORD code = exception_info->ExceptionRecord->ExceptionCode;
    auto reason = fmt::format("Windows Structured Exception: 0x{:08X}", code);

    write_crash_report(reason, boost::stacktrace::stacktrace());

    return EXCEPTION_CONTINUE_SEARCH;
}
#endif

} // anonymous namespace

void install_crash_handlers()
{
    std::set_terminate(terminate_handler);

    std::signal(SIGSEGV, signal_handler);
    std::signal(SIGABRT, signal_handler);
    std::signal(SIGFPE, signal_handler);
    std::signal(SIGILL, signal_handler);
#ifndef _WIN32
    std::signal(SIGBUS, signal_handler);
#endif

#ifdef _WIN32
    SetUnhandledExceptionFilter(unhandled_exception_filter);
#endif
}

} // namespace schrodinger

#endif // __EMSCRIPTEN__
