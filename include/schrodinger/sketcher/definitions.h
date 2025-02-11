// @copyright Schrodinger, LLC - All Rights Reserved

#pragma once

#ifdef _WIN32
#define API_HELPER_IMPORT __declspec(dllimport)
#define API_HELPER_EXPORT __declspec(dllexport)
#else
#define API_HELPER_IMPORT __attribute__((visibility("default")))
#define API_HELPER_EXPORT __attribute__((visibility("default")))
#endif

#ifdef SKETCHER_STATIC_DEFINE
#define SKETCHER_API
#else
#ifdef IN_SKETCHER_DLL
#define SKETCHER_API API_HELPER_EXPORT
#else
#define SKETCHER_API API_HELPER_IMPORT
#endif
#endif
