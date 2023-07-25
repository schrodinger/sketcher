cmake_minimum_required(VERSION 3.16)

project(schrodinger_sketcher VERSION 2023.4 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Library Dependencies
# find_package(Boost REQUIRED COMPONENTS filesystem iostreams serialization system)
find_package(Qt6 REQUIRED COMPONENTS Widgets Svg)
# find_package(coordgen REQUIRED)
# find_package(fmt REQUIRED)
# find_package(maeparser REQUIRED)
# find_package(rdkit REQUIRED)

# TODO: For Qt6.3, use qt_standard_project_setup()
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)
set(CMAKE_AUTOUIC ON)

# Sketcher Library
set(SKETCHER_SRC "${PROJECT_SOURCE_DIR}/src/schrodinger/sketcher")
set(CMAKE_AUTOUIC_SEARCH_PATHS  "${SKETCHER_SRC}/ui")
file(GLOB_RECURSE SOURCE_LIST CONFIGURE_DEPENDS
    ${SKETCHER_SRC}/*.h ${SKETCHER_SRC}/*.cpp)

add_library(schrodinger_sketcher SHARED ${SOURCE_LIST})
target_include_directories(schrodinger_sketcher PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${Qt6Widgets_INCLUDE_DIRS} ${Qt6Svg_INCLUDE_DIRS})
target_link_libraries(schrodinger_sketcher PRIVATE Qt6::Widgets)

# Sketcher Executable
# add_executable(sketcher_app src/bin/app/main.cpp)
# target_include_directories(sketcher_app PRIVATE
#     schrodinger_sketcher ${Qt6Widgets_INCLUDE_DIRS})
# target_link_libraries(sketcher_app PUBLIC
#     schrodinger_sketcher Qt6::Svg Qt6::Widgets)

# Testing Dependencies
find_package(Python)

enable_testing()
add_test (NAME sketcher_tests
  COMMAND ${PYTHON_EXECUTABLE} -m pytest ${CMAKE_CURRENT_SOURCE_DIR}
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/test/
)

################################################################################
# cmake . -B build/ -DCMAKE_PREFIX_PATH=$SCHRODINGER_LIB/Darwin-x86_64/qt-6.2.7/ && cmake --build build/ --clean-first && ctest
################################################################################