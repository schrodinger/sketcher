# Inject <arm_acle.h> into qtbase's qyieldcpu.h.
#
# Apple clang >= 21 reports __has_builtin(__yield) == true, but __yield() is
# only declared in <arm_acle.h>. Without that include the generic branch of
# qYieldCpu() implicitly-declares the function and the Qt build fails under
# -Werror=implicit-function-declaration. Run as a PATCH_COMMAND from the Qt
# source root; idempotent so re-running is harmless.

set(_file "qtbase/src/corelib/thread/qyieldcpu.h")

if(NOT EXISTS "${_file}")
  message(STATUS "qt_yield_acle_patch: ${_file} not found, skipping")
  return()
endif()

file(READ "${_file}" _contents)

if(_contents MATCHES "arm_acle.h")
  message(STATUS "qt_yield_acle_patch: already applied")
  return()
endif()

string(
  REPLACE
    "#include <QtCore/qtconfigmacros.h>"
    "#include <QtCore/qtconfigmacros.h>\n\n#if defined(__has_include)\n#  if defined(__ARM_ACLE) && __has_include(<arm_acle.h>)\n#    include <arm_acle.h>\n#  endif\n#endif"
    _contents
    "${_contents}")

file(WRITE "${_file}" "${_contents}")
message(STATUS "qt_yield_acle_patch: applied to ${_file}")
