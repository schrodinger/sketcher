cmake_minimum_required(VERSION 3.24)

project(
  schrodinger_sketcher
  VERSION 2024.2.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Dependencies
file(READ external/versions.json VERSIONS)
string(JSON BOOST_VERSION GET ${VERSIONS} boost)
string(JSON FMT_VERSION GET ${VERSIONS} fmt)
string(JSON JSON_VERSION GET ${VERSIONS} jsoncpp)
string(JSON QT_VERSION GET ${VERSIONS} qt)
string(JSON RDKIT_VERSION GET ${VERSIONS} rdkit)
string(JSON ZSTD_VERSION GET ${VERSIONS} zstd)

set(Boost_USE_STATIC_LIBS ON) # Ensure static libs are found

find_package(Boost ${BOOST_VERSION} REQUIRED
             COMPONENTS filesystem serialization unit_test_framework)
find_package(Qt6 ${QT_VERSION} REQUIRED COMPONENTS Widgets Svg)
find_package(RDKit ${RDKIT_VERSION} REQUIRED)
find_package(fmt ${FMT_VERSION} REQUIRED)
find_package(jsoncpp REQUIRED) # version not implemented
find_package(zstd ${ZSTD_VERSION} REQUIRED)

# Common library/executable configuration
function(setup_target TARGET)
  if(MSVC)
    target_compile_options(${TARGET} PRIVATE -WX)
  else()
    target_compile_options(${TARGET} PRIVATE -Wall -Werror -Wswitch -O2)
  endif()
  # Allow RDKit includes to be prefixed
  target_include_directories(${TARGET} PRIVATE ${RDKit_INCLUDE_DIRS}/..
                                               ${PROJECT_SOURCE_DIR}/include)
  # Assume use of static libraries
  target_compile_definitions(
    ${TARGET} PRIVATE RDKIT_EXTENSIONS_STATIC_DEFINE SKETCHER_STATIC_DEFINE
                      RDK_BUILD_MAEPARSER_SUPPORT)
endfunction()

# Build static libraries
function(build_static_library TARGET)
  set(INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include/schrodinger/${TARGET})
  set(SRC_DIR ${PROJECT_SOURCE_DIR}/src/schrodinger/${TARGET})
  # Collect public headers so that Qt can build .moc files
  file(GLOB_RECURSE SOURCE_LIST CONFIGURE_DEPENDS ${INCLUDE_DIR}/*.h
       ${SRC_DIR}/*.h ${SRC_DIR}/*.cpp)
  add_library(${TARGET} STATIC ${SOURCE_LIST})
  add_library(schrodinger::${TARGET} ALIAS ${TARGET})
  setup_target(${TARGET})
endfunction()

set(RDKIT_EXTENSIONS_LIB rdkit_extensions)
build_static_library(${RDKIT_EXTENSIONS_LIB})
target_link_libraries(
  ${RDKIT_EXTENSIONS_LIB}
  PRIVATE Boost::filesystem
          Boost::serialization
          RDKit::CIPLabeler
          RDKit::ChemReactions
          RDKit::DetermineBonds
          RDKit::DistGeomHelpers
          RDKit::RDInchiLib
          fmt::fmt
          jsoncpp_lib_static
          zstd::libzstd_static)

set(SKETCHER_LIB sketcher)
build_static_library(${SKETCHER_LIB})
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
set(TEST_DEPENDENCIES Boost::unit_test_framework fmt::fmt Qt6::Widgets
                      ${SKETCHER_LIB})
foreach(TEST_SOURCE ${TEST_SOURCE_LIST})
  cmake_path(GET TEST_SOURCE STEM TEST_TARGET)
  add_executable(${TEST_TARGET} ${TEST_SOURCE})
  setup_target(${TEST_TARGET})
  target_link_libraries(${TEST_TARGET} PRIVATE ${TEST_DEPENDENCIES})
  add_test(NAME ${TEST_TARGET} COMMAND ${TEST_TARGET})
  # Mirror the source directory for test collection
  cmake_path(RELATIVE_PATH TEST_SOURCE OUTPUT_VARIABLE TEST_PATH)
  cmake_path(GET TEST_PATH PARENT_PATH PREFIX_PATH)
  set_target_properties(${TEST_TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                  ${PREFIX_PATH})
endforeach()
