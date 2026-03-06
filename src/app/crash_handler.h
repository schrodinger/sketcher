/* -------------------------------------------------------------------------
 * Crash handler for the standalone schrodinger_sketcher application.
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

#ifndef __EMSCRIPTEN__

namespace schrodinger
{

/**
 * Install crash handlers for unhandled exceptions and fatal signals. Should be
 * called as early as possible in main(), before QApplication. Writes crash
 * reports to a timestamped file in the system temp directory and to stderr.
 */
void install_crash_handlers();

} // namespace schrodinger

#endif // __EMSCRIPTEN__
