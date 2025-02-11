// @copyright Schrodinger, LLC - All Rights Reserved

#pragma once

#define API_HELPER_IMPORT __attribute__((visibility("default")))
#define API_HELPER_EXPORT __attribute__((visibility("default")))

#ifdef SKETCHER_STATIC_DEFINE
#define SKETCHER_API
#else
#ifdef IN_SKETCHER_DLL
#define SKETCHER_API API_HELPER_EXPORT
#else
#define SKETCHER_API API_HELPER_IMPORT
#endif
#endif
