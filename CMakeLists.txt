cmake_minimum_required(VERSION 3.30)

project(
  schrodinger_sketcher
  VERSION 2024.2.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
cmake_policy(SET CMP0167 NEW)

option(BUILD_SHARED_LIBS "Build using shared libraries" OFF)
option(ENFORCE_DEPENDENCY_VERSIONS
    "Require that the dependencies match the versions listed in \
    external/versions.json"
    ON)

set(SUPPRESSION_FILE "test_valgrind.suppressions")
set(MEMTEST_COMMAND
    "valgrind --tool=memcheck --time-stamp=yes --num-callers=20 \
    --gen-suppressions=all --leak-check=yes --keep-debuginfo=no \
    --error-exitcode=29")
separate_arguments(MEMTEST_COMMAND)

if (${ENFORCE_DEPENDENCY_VERSIONS})
    file(READ external/versions.json VERSIONS)
    string(JSON BOOST_VERSION GET ${VERSIONS} boost)
    string(JSON EIGEN_VERSION GET ${VERSIONS} eigen)
    string(JSON FMT_VERSION GET ${VERSIONS} fmt)
    string(JSON QT_VERSION GET ${VERSIONS} qt)
    string(JSON RDKIT_VERSION GET ${VERSIONS} rdkit)
    string(JSON SQLITE_VERSION GET ${VERSIONS} sqlite)
    string(JSON ZLIB_VERSION GET ${VERSIONS} zlib)
    string(JSON ZSTD_VERSION GET ${VERSIONS} zstd)
endif()

if (${BUILD_SHARED_LIBS})
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
find_package(Qt6 ${QT_VERSION} REQUIRED COMPONENTS Widgets Svg)
find_package(Eigen3 ${EIGEN_VERSION} REQUIRED)
find_package(RDKit ${RDKIT_VERSION} REQUIRED)
find_package(fmt ${FMT_VERSION} REQUIRED)
find_package(SQLite3 ${SQLITE_VERSION} REQUIRED CONFIG)
find_package(ZLIB ${ZLIB_VERSION} REQUIRED)
find_package(zstd ${ZSTD_VERSION} REQUIRED)

# Common library/executable configuration
function(setup_target TARGET)
  if(MSVC)
    target_compile_options(${TARGET} PRIVATE -WX)
  elseif(WIN32)
    target_compile_options(${TARGET} PRIVATE -Wall -Werror -Wswitch -O2 -fvisibility=hidden -Wno-restrict)
  else()
    target_compile_options(${TARGET} PRIVATE -Wall -Werror -Wswitch -O2 -fvisibility=hidden)
  endif()
  # Allow RDKit includes to be prefixed
  target_include_directories(${TARGET} PRIVATE ${RDKit_INCLUDE_DIRS}
                                               ${RDKit_INCLUDE_DIRS}/..
                                               ${PROJECT_SOURCE_DIR}/include
                                               ${PROJECT_SOURCE_DIR}/src)
  target_compile_definitions(${TARGET} PRIVATE RDK_BUILD_MAEPARSER_SUPPORT)
  if (NOT ${BUILD_SHARED_LIBS})
    # force use of static libraries
    target_compile_definitions(
        ${TARGET} PRIVATE RDKIT_EXTENSIONS_STATIC_DEFINE SKETCHER_STATIC_DEFINE
                          BOOST_BEAST_HEADER_ONLY)
    if (WIN32)
        target_link_libraries(${TARGET} PUBLIC -static)
    endif()
  endif()
  if (EMSCRIPTEN)
    target_link_options(${TARGET} PRIVATE -sGL_ENABLE_GET_PROC_ADDRESS)
  endif()
endfunction()

# Configure the SKETCHER_LIB and RDKIT_EXTENSIONS_LIB libraries
function(build_schrodinger_library TARGET)
  set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/schrodinger/${TARGET})
  set(SRC_DIR ${PROJECT_SOURCE_DIR}/src/schrodinger/${TARGET})
  # Collect public headers so that Qt can build .moc files
  file(GLOB_RECURSE SOURCE_LIST CONFIGURE_DEPENDS ${INCLUDE_DIR}/*.h
       ${SRC_DIR}/*.h ${SRC_DIR}/*.cpp)
  add_library(${TARGET} ${SOURCE_LIST})
  add_library(schrodinger::${TARGET} ALIAS ${TARGET})
  setup_target(${TARGET})
endfunction()

set(RDKIT_EXTENSIONS_LIB rdkit_extensions)
build_schrodinger_library(${RDKIT_EXTENSIONS_LIB})
target_compile_definitions(${RDKIT_EXTENSIONS_LIB} PRIVATE IN_RDKIT_EXTENSIONS_DLL)
target_link_libraries(
  ${RDKIT_EXTENSIONS_LIB}
  PRIVATE Boost::filesystem
          Boost::iostreams
          Boost::json
          Boost::serialization
          RDKit::CIPLabeler
          RDKit::ChemReactions
          RDKit::DetermineBonds
          RDKit::DistGeomHelpers
          RDKit::MarvinParser
          RDKit::RDInchiLib
          fmt::fmt
          SQLite::SQLite3
          ZLIB::ZLIB
          ${ZSTD_LIB_NAME})

set(SKETCHER_LIB sketcher)
build_schrodinger_library(${SKETCHER_LIB})
target_compile_definitions(${SKETCHER_LIB} PRIVATE IN_SKETCHER_DLL)
target_link_libraries(${SKETCHER_LIB} PRIVATE Qt6::Widgets Qt6::Svg
                                              ${RDKIT_EXTENSIONS_LIB})
# Enable Qt MOC and UIC for the sketcher library
set_target_properties(${SKETCHER_LIB} PROPERTIES AUTOMOC ON)
set_target_properties(${SKETCHER_LIB} PROPERTIES AUTOUIC ON)
set(UI_DIR ${PROJECT_SOURCE_DIR}/src/schrodinger/${SKETCHER_LIB}/ui)
set_target_properties(${SKETCHER_LIB} PROPERTIES AUTOUIC_SEARCH_PATHS ${UI_DIR})

# Standalone executable
set(EXE_TARGET sketcher_app)
add_executable(${EXE_TARGET} src/app/main.cpp)
setup_target(${EXE_TARGET})
target_link_libraries(${EXE_TARGET} PRIVATE Qt6::Widgets ${SKETCHER_LIB})

# Tests
enable_testing()

set(SKETCHER_TEST ${PROJECT_SOURCE_DIR}/test/schrodinger)
file(GLOB_RECURSE TEST_SOURCE_LIST CONFIGURE_DEPENDS ${SKETCHER_TEST}/*.cpp)
set(TEST_DEPENDENCIES Boost::unit_test_framework Eigen3::Eigen fmt::fmt Qt6::Widgets
                      ${RDKIT_EXTENSIONS_LIB} ${SKETCHER_LIB})
file(GLOB GLOBAL_SUPPRESSIONS "test/valgrind_*suppressions*")

foreach(TEST_SOURCE ${TEST_SOURCE_LIST})
  cmake_path(GET TEST_SOURCE STEM TEST_TARGET)
  add_executable(${TEST_TARGET} ${TEST_SOURCE})
  setup_target(${TEST_TARGET})
  target_include_directories(${TEST_TARGET} PRIVATE ${PROJECT_SOURCE_DIR}/test)
  target_link_libraries(${TEST_TARGET} PRIVATE ${TEST_DEPENDENCIES})
  if(EMSCRIPTEN)
    set(TEST_COMMAND node $<TARGET_FILE:${TEST_TARGET}>)
  else()
    set(TEST_COMMAND ${TEST_TARGET})
  endif()
  add_test(NAME ${TEST_TARGET} COMMAND ${TEST_COMMAND})
  # Mirror the source directory for test collection
  cmake_path(RELATIVE_PATH TEST_SOURCE OUTPUT_VARIABLE TEST_PATH)
  cmake_path(GET TEST_PATH PARENT_PATH PREFIX_PATH)
  set_target_properties(${TEST_TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                  ${PREFIX_PATH})
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
    list(APPEND CURRENT_MEMTEST_COMMAND
        "${PREFIX_PATH}/${TEST_COMMAND}")
    add_test(NAME memtest_${TEST_TARGET} COMMAND ${CURRENT_MEMTEST_COMMAND})
    set_property(TEST memtest_${TEST_TARGET} PROPERTY LABELS memtest)
  endif()
endforeach()
