cmake_minimum_required(VERSION 3.16)

project(
  schrodinger_sketcher
  VERSION 2023.4.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_COMPILER clang++) # FIXME: Linux and Windows?
set(CMAKE_OSX_ARCHITECTURES "x86_64")

# FIXME: -WX on Windows
add_definitions(-Wall -Werror -Wswitch -O2 -fvisibility=hidden)

# Library Dependencies
find_package(Qt6 6.2.7 REQUIRED COMPONENTS Widgets Svg)
find_package(rdkit 2023.03.2 REQUIRED)
find_package(maeparser REQUIRED) # FIXME: add version support SHARED-10008
find_package(coordgen REQUIRED)
find_package(Boost 1.76.0 REQUIRED COMPONENTS system unit_test_framework)
find_package(fmt 6.1.2 REQUIRED)

set(SKETCHER_DEP_INCLUDES
    ${Qt6Widgets_INCLUDE_DIRS} ${Qt6Svg_INCLUDE_DIRS}
    ${rdkit_ROOT}/include/rdkit ${coordgen_ROOT}/include ${Boost_INCLUDE_DIRS}
    ${FMT_ROOT}/lib/cmake/include)
set(SKETCHER_DEP_LIBS
    Qt6::Widgets
    Qt6::Svg
    RDKit::CIPLabeler
    RDKit::RDInchiLib
    RDKit::ChemReactions
    Boost::system
    fmt::fmt)

# TODO: For Qt6.3, use qt_standard_project_setup()
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

set(SKETCHER_INC ${PROJECT_SOURCE_DIR}/include/schrodinger/sketcher)
set(SKETCHER_SRC ${PROJECT_SOURCE_DIR}/src/schrodinger/sketcher)
set(SKETCHER_TEST ${PROJECT_SOURCE_DIR}/test/schrodinger/sketcher)

# Sketcher Shared Library
set(LIB_NAME ${CMAKE_PROJECT_NAME})
file(GLOB_RECURSE SOURCE_LIST CONFIGURE_DEPENDS ${SKETCHER_INC}/*.h
     ${SKETCHER_SRC}/*.h ${SKETCHER_SRC}/*.cpp)
set(CMAKE_AUTOUIC_SEARCH_PATHS ${SKETCHER_SRC}/ui)
set(SKETCHER_LIB_INCLUDES ${PROJECT_SOURCE_DIR}/include
                          ${SKETCHER_DEP_INCLUDES})
add_library(${LIB_NAME} SHARED ${SOURCE_LIST})
set_target_properties(${LIB_NAME} PROPERTIES DEFINE_SYMBOL IN_SKETCHER_DLL)
target_include_directories(${LIB_NAME} PRIVATE ${SKETCHER_LIB_INCLUDES})
target_link_libraries(${LIB_NAME} PRIVATE ${SKETCHER_DEP_LIBS})
target_compile_definitions(${LIB_NAME} PRIVATE "-DRDK_BUILD_MAEPARSER_SUPPORT")

# Sketcher Executable
set(EXE_NAME sketcher_app)
add_executable(${EXE_NAME} src/bin/app/main.cpp)
set(EXE_DEP_INCLUDES ${PROJECT_SOURCE_DIR}/include ${Qt6Widgets_INCLUDE_DIRS})
target_include_directories(${EXE_NAME} PRIVATE ${EXE_DEP_INCLUDES})
target_link_libraries(${EXE_NAME} PRIVATE ${LIB_NAME} Qt6::Widgets)

# TODO: Sketcher WASM

# Tests
enable_testing()

file(GLOB_RECURSE TEST_SOURCE_LIST CONFIGURE_DEPENDS ${SKETCHER_TEST}/*.cpp)
set(TEST_DEP_LIBS ${LIB_NAME} Boost::unit_test_framework ${SKETCHER_DEP_LIBS})
foreach(TEST_SOURCE ${TEST_SOURCE_LIST})
  get_filename_component(TEST_NAME ${TEST_SOURCE} NAME_WE)
  add_executable(${TEST_NAME} ${TEST_SOURCE})
  target_include_directories(${TEST_NAME} PRIVATE ${SKETCHER_LIB_INCLUDES})
  target_link_libraries(${TEST_NAME} PRIVATE ${TEST_DEP_LIBS})
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
endforeach()
