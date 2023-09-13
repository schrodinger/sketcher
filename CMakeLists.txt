cmake_minimum_required(VERSION 3.20)

project(
  schrodinger_sketcher
  VERSION 2023.4.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# WASM specific settings
if(CMAKE_CXX_COMPILER MATCHES "/em\\+\\+(-[a-zA-Z0-9.])?$")
  set(WASM_BUILD ON)
  set(CMAKE_CXX_COMPILER_ID "Emscripten")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s DISABLE_EXCEPTION_CATCHING=0")
endif()

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}")

# FIXME: -WX on Windows
add_definitions(-Wall -Werror -Wswitch -O2 -fvisibility=hidden)

# Library Dependencies
find_package(Boost 1.76.0 REQUIRED COMPONENTS iostreams serialization system
                                              unit_test_framework)
find_package(Qt6 6.2.7 REQUIRED COMPONENTS Widgets Svg)
find_package(fmt 6.1.2 REQUIRED)
find_package(maeparser REQUIRED) # FIXME: add version support SHARED-10008
find_package(rdkit 2023.03.2 REQUIRED)

if(NOT ${RDKit_HAS_MAEPARSER_SUPPORT} STREQUAL "ON")
  message(FATAL_ERROR "RDKit must be built with maeparser support")
endif()

set(SKETCHER_DEP_INCLUDES
    ${Boost_INCLUDE_DIRS} ${FMT_ROOT}/lib/cmake/include ${Qt6Svg_INCLUDE_DIRS}
    ${Qt6Widgets_INCLUDE_DIRS} ${rdkit_ROOT}/include/rdkit)
set(SKETCHER_DEP_LIBS
    Boost::iostreams
    Boost::serialization
    Boost::system
    Qt6::Svg
    Qt6::Widgets
    RDKit::CIPLabeler
    RDKit::ChemReactions
    RDKit::RDInchiLib
    fmt::fmt)

# TODO: For Qt6.3, use qt_standard_project_setup()
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

set(SKETCHER_INC ${PROJECT_SOURCE_DIR}/include/schrodinger/sketcher)
set(SKETCHER_SRC ${PROJECT_SOURCE_DIR}/src/schrodinger/sketcher)
set(SKETCHER_TEST ${PROJECT_SOURCE_DIR}/test/schrodinger/sketcher)

# Sketcher Shared Library
set(LIB_TARGET ${CMAKE_PROJECT_NAME})
file(GLOB_RECURSE SOURCE_LIST CONFIGURE_DEPENDS ${SKETCHER_INC}/*.h
     ${SKETCHER_SRC}/*.h ${SKETCHER_SRC}/*.cpp)
set(CMAKE_AUTOUIC_SEARCH_PATHS ${SKETCHER_SRC}/ui)
set(SKETCHER_LIB_INCLUDES ${PROJECT_SOURCE_DIR}/include
                          ${SKETCHER_DEP_INCLUDES})
add_library(${LIB_TARGET} SHARED ${SOURCE_LIST})
set_target_properties(${LIB_TARGET} PROPERTIES DEFINE_SYMBOL IN_SKETCHER_DLL)
target_include_directories(${LIB_TARGET} PRIVATE ${SKETCHER_LIB_INCLUDES})
target_link_libraries(${LIB_TARGET} PRIVATE ${SKETCHER_DEP_LIBS})
target_compile_definitions(${LIB_TARGET}
                           PRIVATE "-DRDK_BUILD_MAEPARSER_SUPPORT")

# Sketcher Executable
set(EXE_TARGET sketcher_app)
add_executable(${EXE_TARGET} src/app/main.cpp)
set(EXE_DEP_INCLUDES ${PROJECT_SOURCE_DIR}/include ${Qt6Widgets_INCLUDE_DIRS})
target_include_directories(${EXE_TARGET} PRIVATE ${EXE_DEP_INCLUDES})
target_link_libraries(${EXE_TARGET} PRIVATE ${LIB_TARGET} Qt6::Widgets)

# TODO: Sketcher WASM

# Tests
enable_testing()

file(GLOB_RECURSE TEST_SOURCE_LIST CONFIGURE_DEPENDS ${SKETCHER_TEST}/*.cpp)
set(TEST_DEP_LIBS ${LIB_TARGET} Boost::unit_test_framework ${SKETCHER_DEP_LIBS})
foreach(TEST_SOURCE ${TEST_SOURCE_LIST})
  cmake_path(GET TEST_SOURCE STEM TEST_TARGET)
  add_executable(${TEST_TARGET} ${TEST_SOURCE})
  target_include_directories(${TEST_TARGET} PRIVATE ${SKETCHER_LIB_INCLUDES})
  target_link_libraries(${TEST_TARGET} PRIVATE ${TEST_DEP_LIBS})
  add_test(NAME ${TEST_TARGET} COMMAND ${TEST_TARGET})
  # Mirror the source directory for test collection
  cmake_path(RELATIVE_PATH TEST_SOURCE OUTPUT_VARIABLE TEST_PATH)
  cmake_path(GET TEST_PATH PARENT_PATH PREFIX_PATH)
  set_target_properties(${TEST_TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                  ${PREFIX_PATH})
endforeach()
