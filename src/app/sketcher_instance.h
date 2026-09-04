/* -------------------------------------------------------------------------
 * Access to the application's single SketcherWidget.
 *
 * Defined in main.cpp; declared here so that other translation units in the
 * app can reach the instance without duplicating the declaration.
 *
 * Copyright Schrodinger LLC, All Rights Reserved.
 --------------------------------------------------------------------------- */

#pragma once

namespace schrodinger
{
namespace sketcher
{
class SketcherWidget;
} // namespace sketcher
} // namespace schrodinger

schrodinger::sketcher::SketcherWidget& get_sketcher_instance();
