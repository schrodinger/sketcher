// @copyright Schrodinger, LLC - All Rights Reserved

#pragma once

#define API_HELPER_IMPORT __attribute__((visibility("default")))
#define API_HELPER_EXPORT __attribute__((visibility("default")))

#ifdef RDKIT_EXTENSIONS_STATIC_DEFINE
#define RDKIT_EXTENSIONS_API
#else
#ifdef IN_RDKIT_EXTENSIONS_DLL
#define RDKIT_EXTENSIONS_API API_HELPER_EXPORT
#else
#define RDKIT_EXTENSIONS_API API_HELPER_IMPORT
#endif
#endif
