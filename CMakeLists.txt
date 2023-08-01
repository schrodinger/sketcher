cmake_minimum_required(VERSION 3.20)

project(
  schrodinger_sketcher
  VERSION 2023.4.0
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID}")

# FIXME: -WX on Windows
add_definitions(-Wall -Werror -Wswitch -O2 -fvisibility=hidden)

# Sketcher Executable
set(EXE_TARGET sketcher_app)
add_executable(${EXE_TARGET} src/app/main.cpp)
