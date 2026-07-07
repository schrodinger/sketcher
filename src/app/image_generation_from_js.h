/**
 * Helper functions for calling the image_generation API from JavaScript
 */

#pragma once

#ifdef __EMSCRIPTEN__

#include <emscripten/val.h>

#include <QByteArray>

#include "schrodinger/sketcher/image_generation.h"

/**
 * Convert a QByteArray to a JavaScript array of uint8 values
 */
emscripten::val qbyte_array_to_uint8_array(QByteArray bytes);

/**
 * Convert a JavaScript object that describes render options to an actual C++
 * RenderOptions instance.
 */
schrodinger::sketcher::RenderOptions
render_options_from_js(const emscripten::val& options);

#endif
