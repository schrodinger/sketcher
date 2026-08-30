// @copyright Schrodinger, LLC - All Rights Reserved

#pragma once

#ifndef API_HELPER_IMPORT
#ifdef _WIN32
#define API_HELPER_IMPORT __declspec(dllimport)
#else
#define API_HELPER_IMPORT __attribute__((visibility("default")))
#endif
#endif

#ifndef API_HELPER_EXPORT
#ifdef _WIN32
#define API_HELPER_EXPORT __declspec(dllexport)
#else
#define API_HELPER_EXPORT __attribute__((visibility("default")))
#endif
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
