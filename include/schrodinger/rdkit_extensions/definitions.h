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

#ifdef RDKIT_EXTENSIONS_STATIC_DEFINE
#define RDKIT_EXTENSIONS_API
#else
#ifdef IN_RDKIT_EXTENSIONS_DLL
#define RDKIT_EXTENSIONS_API API_HELPER_EXPORT
#else
#define RDKIT_EXTENSIONS_API API_HELPER_IMPORT
#endif
#endif
