cmake_minimum_required(VERSION 3.16)

project(schrodinger_sketcher VERSION 2023.4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_COMPILER clang++) # FIXME: Linux and Windows?
set(CMAKE_OSX_ARCHITECTURES "x86_64")

# Library Dependencies
find_package(Qt6 REQUIRED COMPONENTS Widgets Svg)
# find_package(rdkit REQUIRED)
# find_package(coordgen REQUIRED)
# find_package(maeparser REQUIRED)
find_package(Boost REQUIRED COMPONENTS
    filesystem iostreams serialization system unit_test_framework)
find_package(fmt REQUIRED)

# TODO: For Qt6.3, use qt_standard_project_setup()
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# FIXME: -WX on Windows
add_definitions(-Wall -Werror -Wswitch -O2 -fvisibility=hidden)

# Sketcher Library
set(SKETCHER_SRC "${PROJECT_SOURCE_DIR}/src/schrodinger/sketcher")
set(CMAKE_AUTOUIC_SEARCH_PATHS  "${SKETCHER_SRC}/ui")
file(GLOB_RECURSE SOURCE_LIST CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/include/schrodinger/sketcher/*.h
    ${SKETCHER_SRC}/*.h ${SKETCHER_SRC}/*.cpp)

add_library(schrodinger_sketcher SHARED ${SOURCE_LIST})
set_target_properties(schrodinger_sketcher PROPERTIES
    DEFINE_SYMBOL IN_SKETCHER_DLL)
target_include_directories(schrodinger_sketcher PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${Qt6Widgets_INCLUDE_DIRS} ${Qt6Svg_INCLUDE_DIRS}
    # rdkit coordgen maeparser
    ${Boost_INCLUDE_DIRS}
    ${FMT_ROOT}/lib/cmake/include)
target_link_libraries(schrodinger_sketcher PRIVATE
    Qt6::Widgets Qt6::Svg Boost::filesystem fmt::fmt)

# Sketcher Executable
add_executable(sketcher_app src/bin/app/main.cpp)
target_include_directories(sketcher_app PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${Qt6Widgets_INCLUDE_DIRS})
target_link_libraries(sketcher_app PRIVATE schrodinger_sketcher Qt6::Widgets)

# TODO: Sketcher WASM

# Tests
enable_testing()

file(GLOB_RECURSE TEST_SOURCE_LIST CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/test/schrodinger/sketcher/*.cpp)
foreach(test_source ${TEST_SOURCE_LIST})
    get_filename_component(test_name ${test_source} NAME_WE)
    add_executable(${test_name} ${test_source})
    target_include_directories(${test_name} PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${Qt6Widgets_INCLUDE_DIRS}
        # rdkit coordgen maeparser
        ${Boost_INCLUDE_DIRS}
        ${FMT_ROOT}/lib/cmake/include)
    target_link_libraries(${test_name} PRIVATE
        schrodinger_sketcher Qt6::Widgets Boost::unit_test_framework fmt::fmt)
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()


################################################################################
# rm -rf build/* && cmake . -B build/ -DCMAKE_PREFIX_PATH="$SCHRODINGER_LIB/Darwin-x86_64/boost_1_81_0-py38_cpp20/;$SCHRODINGER_LIB/Darwin-x86_64/qt-6.2.7/;$SCHRODINGER_LIB/Darwin-x86_64/fmt-6.1.2" && cmake --build build/ && ctest --test-dir build/
################################################################################