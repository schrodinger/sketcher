cmake_minimum_required(VERSION 3.24)

file(READ ${CMAKE_CURRENT_SOURCE_DIR}/version.json VERSIONS)
string(JSON SKETCHER_VERSION GET ${VERSIONS} sketcher)

project(
  schrodinger_sketcher
  VERSION ${SKETCHER_VERSION}
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
option(ENABLE_TESTING "Build and enable tests" ON)
option(ENFORCE_DEPENDENCY_VERSIONS "Require versions.json enforcement" ON)

# Load dependency versions
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/external/versions.json VERSIONS)
if(${ENFORCE_DEPENDENCY_VERSIONS})
  string(JSON BOOST_VERSION GET ${VERSIONS} boost)
  string(JSON EIGEN_VERSION GET ${VERSIONS} eigen)
  string(JSON FMT_VERSION GET ${VERSIONS} fmt)
  string(JSON QT_VERSION GET ${VERSIONS} qt)
  string(JSON RDKIT_VERSION GET ${VERSIONS} rdkit)
  string(JSON SQLITE_VERSION GET ${VERSIONS} sqlite)
  string(JSON ZLIB_VERSION GET ${VERSIONS} zlib)
  string(JSON ZSTD_VERSION GET ${VERSIONS} zstd)
endif()

if(${BUILD_SHARED_LIBS})
  set(ZSTD_LIB_NAME zstd::libzstd_shared)
  set(Boost_USE_STATIC_LIBS OFF)
  # RDKit needs the shared versions of maeparser and coordgen
  find_package(maeparser REQUIRED)
  find_package(coordgen REQUIRED)
else()
  # we're building a static build
  set(ZSTD_LIB_NAME zstd::libzstd_static)
  set(Boost_USE_STATIC_LIBS ON)
endif()

if(CMAKE_FIND_ROOT_PATH)
  # make sure we can always find the dependencies that we've built, even if
  # emscripten is trying to prevent us from finding system libraries
  list(APPEND CMAKE_FIND_ROOT_PATH ${CMAKE_CURRENT_BINARY_DIR})
endif()

if(EMSCRIPTEN)
  set(QT_HOST_PATH ${CMAKE_CURRENT_BINARY_DIR}/external/host/qt-${QT_VERSION})
endif()

# we don't include boost::beast as a required component since GitHub builds use
# it as a non-built, headers-only library, which means that there's no
# associated config.cmake file for find_package to recognize
find_package(
  Boost ${BOOST_VERSION} REQUIRED COMPONENTS filesystem iostreams json
                                             serialization unit_test_framework)
find_package(Qt6 ${QT_VERSION} REQUIRED COMPONENTS Widgets Svg SvgWidgets Test)
find_package(Eigen3 ${EIGEN_VERSION} REQUIRED)
find_package(RDKit ${RDKIT_VERSION} REQUIRED)
find_package(fmt ${FMT_VERSION} REQUIRED)
find_package(SQLite3 ${SQLITE_VERSION} REQUIRED)
find_package(ZLIB ${ZLIB_VERSION} REQUIRED)
find_package(zstd ${ZSTD_VERSION} REQUIRED)

# Common library/executable configuration
function(setup_target TARGET)
  # Visibility settings for dynamic vs static builds
  if(${BUILD_SHARED_LIBS})
    set(VISIBILITY "default")
  else()
    set(VISIBILITY "hidden")
  endif()
  # Compilation optimization, warnings, and error handling
  if(MSVC)
    target_compile_options(${TARGET} PRIVATE -WX -EHsc)
    if(${BUILD_SHARED_LIBS})
      target_compile_definitions(${TARGET} PRIVATE RDKIT_DYN_LINK)
    endif()
  elseif(WIN32)
    target_compile_options(
      ${TARGET} PRIVATE -Wall -Werror -Wswitch -O2 -fvisibility=${VISIBILITY}
                        -Wno-restrict)
  else()
    target_compile_options(${TARGET} PRIVATE -Wall -Werror -Wswitch -O2
                                             -fvisibility=${VISIBILITY})
  endif()
  # Allow RDKit includes to be prefixed, and add all relevant and necessary
  # sketcher-specific include directories
  target_include_directories(
    ${TARGET}
    PRIVATE ${RDKit_INCLUDE_DIRS} ${RDKit_INCLUDE_DIRS}/..
            ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include
            ${PROJECT_SOURCE_DIR}/src)
  target_compile_definitions(${TARGET} PRIVATE RDK_BUILD_MAEPARSER_SUPPORT)
  if(NOT ${BUILD_SHARED_LIBS})
    # force use of static libraries
    target_compile_definitions(
      ${TARGET} PRIVATE RDKIT_EXTENSIONS_STATIC_DEFINE SKETCHER_STATIC_DEFINE
                        BOOST_BEAST_HEADER_ONLY)
    if(WIN32)
      if(MSVC)
        # use the static version of the runtime library instead of the default
        # shared version
        set_property(
          TARGET ${TARGET} PROPERTY MSVC_RUNTIME_LIBRARY
                                    "MultiThreaded$<$<CONFIG:Debug>:Debug>")
      else()
        target_link_options(${TARGET} PUBLIC -static)
      endif()
    endif()
  endif()
  if(EMSCRIPTEN)
    # We optimize for space in WASM builds because otherwise binaryen (run
    # because of ASYNCIFY) runs out of memory
    target_link_options(${TARGET} PRIVATE -sGL_ENABLE_GET_PROC_ADDRESS
                        -fexceptions -sDISABLE_EXCEPTION_CATCHING=0 -Os)
    target_compile_options(${TARGET} PRIVATE -fexceptions
                                             -sDISABLE_EXCEPTION_CATCHING=0 -Os)
  endif()
endfunction()

# Configure the SKETCHER_LIB and RDKIT_EXTENSIONS_LIB libraries
function(build_schrodinger_library TARGET)
  set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/schrodinger/${TARGET})
  set(SRC_DIR ${PROJECT_SOURCE_DIR}/src/schrodinger/${TARGET})
  # Collect all headers, source, and resources; note we collect public headers
  # so that Qt can build .moc files
  file(
    GLOB_RECURSE
    SOURCE_LIST
    CONFIGURE_DEPENDS
    ${INCLUDE_DIR}/*.h
    ${SRC_DIR}/*.h
    ${SRC_DIR}/*.cpp
    ${SRC_DIR}/*.qrc)
  add_library(${TARGET} ${SOURCE_LIST})
  add_library(schrodinger::${TARGET} ALIAS ${TARGET})
  set_target_properties(${TARGET} PROPERTIES OUTPUT_NAME schrodinger_${TARGET})
  setup_target(${TARGET})
  # Auto-generate compile definitions
  string(TOUPPER ${TARGET} TARGET_UPPER)
  target_compile_definitions(${TARGET} PRIVATE IN_${TARGET_UPPER}_DLL)
endfunction()

# Enable Qt's MOC and UIC features
function(enable_qt_moc TARGET)
  set(UI_DIR ${PROJECT_SOURCE_DIR}/src/schrodinger/sketcher/ui)
  set_target_properties(
    ${TARGET}
    PROPERTIES AUTOMOC ON
               AUTOUIC ON
               AUTOUIC_SEARCH_PATHS ${UI_DIR}
               AUTORCC ON)
endfunction()

# Add a dependency on the generated header
function(add_generated_header_dependency TARGET)
  # The header is generated with a Python script, so we need to find the Python
  # interpreter
  find_package(
    Python
    COMPONENTS Interpreter
    REQUIRED)
  set(SRC_DIR ${PROJECT_SOURCE_DIR}/src/schrodinger/rdkit_extensions)
  set(HDR_PATH
      ${PROJECT_BINARY_DIR}/include/schrodinger/rdkit_extensions/monomer_db_default_monomers.h
  )
  add_custom_command(
    OUTPUT ${HDR_PATH}
    COMMAND ${Python_EXECUTABLE} ${SRC_DIR}/generate_default_monomers_header.py
            ${SRC_DIR}/default_monomer_definitions.json ${HDR_PATH}
    DEPENDS ${SRC_DIR}/default_monomer_definitions.json
            ${SRC_DIR}/generate_default_monomers_header.py
    COMMENT "Generating monomer_db_default_monomers.h")
  add_custom_target(generated_header DEPENDS ${HDR_PATH})
  add_dependencies(${TARGET} generated_header)
endfunction()

# schrodinger::rdkit_extensions library
set(RDKIT_EXTENSIONS_TARGET rdkit_extensions)
build_schrodinger_library(${RDKIT_EXTENSIONS_TARGET})
add_generated_header_dependency(${RDKIT_EXTENSIONS_TARGET})
target_link_libraries(
  ${RDKIT_EXTENSIONS_TARGET}
  PRIVATE Boost::filesystem
          Boost::iostreams
          Boost::json
          Boost::serialization
          fmt::fmt
          RDKit::ChemReactions
          RDKit::DetermineBonds
          RDKit::DistGeomHelpers
          RDKit::MarvinParser
          RDKit::RDInchiLib
          SQLite::SQLite3
          ZLIB::ZLIB
          ${ZSTD_LIB_NAME})

# schrodinger::sketcher library
set(SKETCHER_TARGET sketcher)
build_schrodinger_library(${SKETCHER_TARGET})
target_link_libraries(
  ${SKETCHER_TARGET}
  PRIVATE Boost::boost
          Boost::filesystem
          fmt::fmt
          Qt6::Svg
          Qt6::SvgWidgets
          Qt6::Widgets
          RDKit::CIPLabeler
          RDKit::ChemReactions
          ${RDKIT_EXTENSIONS_TARGET})
# Configure and include version information for the sketcher
configure_file(${PROJECT_SOURCE_DIR}/include/schrodinger/sketcher/version.h.in
               ${PROJECT_BINARY_DIR}/include/schrodinger/sketcher/version.h)
enable_qt_moc(${SKETCHER_TARGET})

# Standalone executable
set(APP_TARGET schrodinger_sketcher)
add_executable(${APP_TARGET} src/app/main.cpp)
setup_target(${APP_TARGET})
if(EMSCRIPTEN)
  set(EXTRA_SKETCHER_APP_LINK_FLAGS
      -sEXPORT_NAME=schrodinger_sketcher_entry
      -sMODULARIZE=1
      -sEXPORTED_RUNTIME_METHODS=[callMain,FS,JSEvents,specialHTMLTargets,stringToUTF16,UTF16ToString]
      -sASYNCIFY
      -lembind
      -fexceptions
      -sDISABLE_EXCEPTION_CATCHING=0)
endif()
target_link_libraries(
  ${APP_TARGET}
  PRIVATE Boost::boost Boost::filesystem Qt6::Widgets
          ${RDKIT_EXTENSIONS_TARGET} ${SKETCHER_TARGET}
          ${EXTRA_SKETCHER_APP_LINK_FLAGS})
if(MSVC AND NOT ${BUILD_SHARED_LIBS})
  # When using MSVC and building static libs, this target generates a file named
  # schrodinger_sketcher.lib, which conflicts with an identically named file
  # from the SKETCHER_TARGET, leading to an error from the linker. By setting
  # IMPORT_SUFFIX, we change the filename for this target to
  # schrodinger_sketcher_bin.lib.
  set_target_properties(
    ${APP_TARGET} PROPERTIES IMPORT_SUFFIX _bin${CMAKE_IMPORT_LIBRARY_SUFFIX})
endif()

# Tests
if(ENABLE_TESTING)
  enable_testing()

  set(SKETCHER_TEST ${PROJECT_SOURCE_DIR}/test/schrodinger)
  file(GLOB_RECURSE TEST_SOURCE_LIST CONFIGURE_DEPENDS ${SKETCHER_TEST}/*.cpp)
  set(TEST_DEPENDENCIES
      Boost::filesystem
      Boost::unit_test_framework
      Eigen3::Eigen
      fmt::fmt
      Qt6::Test
      Qt6::Widgets
      RDKit::ChemReactions
      RDKit::FileParsers
      ${RDKIT_EXTENSIONS_TARGET}
      ${SKETCHER_TARGET})

  file(GLOB GLOBAL_SUPPRESSIONS "test/valgrind_*suppressions*")
  set(SUPPRESSION_FILE "test_valgrind.suppressions")
  set(MEMTEST_COMMAND
      "valgrind --tool=memcheck --time-stamp=yes --num-callers=20 \
        --gen-suppressions=all --leak-check=full \
        --keep-debuginfo=no --error-exitcode=29")
  separate_arguments(MEMTEST_COMMAND)

  # List of tests that require Qt MOC/UIC features (Q_OBJECT or UI files)
  set(TESTS_REQUIRING_QT_MOC
      test_abstract_undoable_model
      test_atom_query_popup
      test_bond_order_popup
      test_bond_query_popup
      test_bracket_subgroup_dialog
      test_draw_tools_widget
      test_edit_atom_properties
      test_enumeration_tool_widget
      test_file_export_dialog
      test_paste_in_text_dialog
      test_periodic_table_widget
      test_ring_tool_widget
      test_select_options_widget
      test_set_atom_widget
      test_sketcher_top_bar
      test_sketcher_widget
      test_stereo_bond_popup)

  foreach(TEST_SOURCE ${TEST_SOURCE_LIST})
    cmake_path(GET TEST_SOURCE STEM TEST_STEM)
    cmake_path(RELATIVE_PATH TEST_SOURCE OUTPUT_VARIABLE TEST_PATH)
    cmake_path(GET TEST_PATH PARENT_PATH PREFIX_PATH)
    # Create unique target name by replacing path separators
    string(REPLACE "/" "-" TEST_PREFIX "${PREFIX_PATH}")
    set(TEST_TARGET "${TEST_PREFIX}-${TEST_STEM}")
    add_executable(${TEST_TARGET} ${TEST_SOURCE})
    setup_target(${TEST_TARGET})
    target_include_directories(${TEST_TARGET}
                               PRIVATE ${PROJECT_SOURCE_DIR}/test)
    target_link_libraries(${TEST_TARGET} PRIVATE ${TEST_DEPENDENCIES})
    if(MSVC)
      # ignore discarded [[nodiscard]] values for the tests
      target_compile_options(${TEST_TARGET} PRIVATE -wd4834)
    endif()

    if(TEST_STEM IN_LIST TESTS_REQUIRING_QT_MOC)
      enable_qt_moc(${TEST_TARGET})
    endif()

    if(EMSCRIPTEN)
      set(TEST_COMMAND node $<TARGET_FILE:${TEST_TARGET}>)
    else()
      set(TEST_COMMAND ${TEST_TARGET})
    endif()
    set_target_properties(
      ${TEST_TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PREFIX_PATH}
                                OUTPUT_NAME ${TEST_STEM})
    # Create memtest target and command
    if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
      # create the valgrind test
      set(CURRENT_MEMTEST_COMMAND ${MEMTEST_COMMAND})
      if(EXISTS "${PREFIX_PATH}/${SUPPRESSION_FILE}")
        list(APPEND CURRENT_MEMTEST_COMMAND
             "--suppressions=${PREFIX_PATH}/${SUPPRESSION_FILE}")
      endif()
      foreach(CURRENT_GLOBAL_SUPPRESSION ${GLOBAL_SUPPRESSIONS})
        list(APPEND CURRENT_MEMTEST_COMMAND
             "--suppressions=${CURRENT_GLOBAL_SUPPRESSION}")
      endforeach()
      list(APPEND CURRENT_MEMTEST_COMMAND "$<TARGET_FILE:${TEST_TARGET}>")
      set(TEST_TARGET "memtest_${TEST_TARGET}")
      set(TEST_COMMAND "${CURRENT_MEMTEST_COMMAND}")
    endif()
    # Add the test to CTest
    add_test(NAME ${TEST_TARGET} COMMAND ${TEST_COMMAND})
    set_tests_properties(
      ${TEST_TARGET} PROPERTIES ENVIRONMENT
                                "SKETCHER_SOURCE_DIR=${CMAKE_SOURCE_DIR}")
    # Add label for memtest targets
    if(${CMAKE_BUILD_TYPE} STREQUAL "Debug")
      set_property(TEST ${TEST_TARGET} PROPERTY LABELS memtest)
    endif()
  endforeach()
endif()

# Install public headers, generated version header, libraries, and executable
install(DIRECTORY ${PROJECT_SOURCE_DIR}/include/schrodinger
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(FILES ${PROJECT_BINARY_DIR}/include/schrodinger/sketcher/version.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/schrodinger/sketcher)
install(TARGETS ${APP_TARGET} DESTINATION ${CMAKE_INSTALL_BINDIR})
install(TARGETS ${SKETCHER_TARGET} ${RDKIT_EXTENSIONS_TARGET}
        DESTINATION ${CMAKE_INSTALL_LIBDIR})
