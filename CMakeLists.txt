cmake_minimum_required(VERSION 3.20)

# Versions
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/versions.json VERSION_JSON)
string(JSON SKETCHER_VERSION GET ${VERSION_JSON} schrodinger_sketcher)
string(JSON FMT_VERSION GET ${VERSION_JSON} fmt)

project(
  schrodinger_sketcher
  VERSION ${SKETCHER_VERSION}
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}")

# FIXME: -WX on Windows
add_definitions(-Wall -Werror -Wswitch -O2 -fvisibility=hidden)

# Library Dependencies
find_package(fmt ${FMT_VERSION} REQUIRED)

set(SKETCHER_DEP_INCLUDES ${FMT_ROOT}/lib/cmake/include)
set(SKETCHER_DEP_LIBS fmt::fmt)

# Sketcher Executable
set(EXE_TARGET sketcher_app)
add_executable(${EXE_TARGET} src/app/main.cpp)
target_include_directories(${EXE_TARGET} PRIVATE ${SKETCHER_DEP_INCLUDES})
target_link_libraries(${EXE_TARGET} PRIVATE ${SKETCHER_DEP_LIBS})
