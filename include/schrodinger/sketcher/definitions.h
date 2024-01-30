/**
 * Public definitions for schrodinger::sketcher:: library
 */

#pragma once

#ifdef SKETCHER_BUILD_SHARED_LIBS
#ifdef _WIN32
#define API_HELPER_IMPORT __declspec(dllimport)
#define API_HELPER_EXPORT __declspec(dllexport)
#else
#define API_HELPER_IMPORT __attribute__((visibility("default")))
#define API_HELPER_EXPORT __attribute__((visibility("default")))
#endif
#endif

#ifdef IN_SKETCHER_DLL
#define SKETCHER_API API_HELPER_EXPORT
#else
#define SKETCHER_API API_HELPER_IMPORT
#endif
