cmake_minimum_required(VERSION 3.24)
include(FetchContent)

project(
  schrodinger_sketcher
  VERSION 2024.2.0
  LANGUAGES CXX)

set(BOOST_VERSION 1.81.0)
set(FMT_VERSION 10.1.0)
set(ZSTD_VERSION 1.5.5)

# C++ Configurartion
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Boost Dependency set(Boost_DEBUG ON)
set(BOOST_ENABLE_CMAKE ON)
set(Boost_USE_SHARED_LIBS OFF)
set(BOOST_INCLUDE_LIBRARIES
    assign
    beast
    bimap
    filesystem
    serialization
    system
    test)
set(BOOST_PACKAGE_COMPONENTS filesystem serialization system
                             unit_test_framework)
FetchContent_Declare(
  Boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG boost-${BOOST_VERSION}
  GIT_SHALLOW TRUE
  FIND_PACKAGE_ARGS ${BOOST_VERSION} COMPONENTS ${BOOST_PACKAGE_COMPONENTS})

# fmt Dependency
set(FMT_DOC OFF)
set(FMT_TEST OFF)
set(FMT_INSTALL OFF)
FetchContent_Declare(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG ${FMT_VERSION}
  GIT_SHALLOW TRUE
  FIND_PACKAGE_ARGS ${FMT_VERSION})

# zstd Dependency
set(ZSTD_BUILD_SHARED OFF)
set(ZSTD_BUILD_PROGRAMS OFF)
FetchContent_Declare(
  zstd
  GIT_REPOSITORY https://github.com/facebook/zstd.git
  GIT_TAG v${ZSTD_VERSION}
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR build/cmake FIND_PACKAGE_ARGS ${ZSTD_VERSION})

# Fetch all missing dependencies before build
FetchContent_MakeAvailable(Boost fmt zstd)

# Platform-specific compiler flags
if(MSVC)
  set(SKETCHER_COMPILE_FLAGS -WX)
else()
  set(SKETCHER_COMPILE_FLAGS -Wall -Werror -Wswitch -fvisibility=hidden)
endif()

# Source paths
set(SKETCHER_INC ${PROJECT_SOURCE_DIR}/include/schrodinger/sketcher)
set(SKETCHER_SRC ${PROJECT_SOURCE_DIR}/src/schrodinger/sketcher)
set(SKETCHER_TEST ${PROJECT_SOURCE_DIR}/test/schrodinger/sketcher)

set(SKETCHER_LIB_INCLUDES
    ${PROJECT_SOURCE_DIR}/include
    # zstd target lacks includes https://github.com/facebook/zstd/issues/3758
    ${zstd_SOURCE_DIR}/lib)

set(SKETCHER_DEPS Boost::filesystem Boost::serialization Boost::system fmt::fmt
                  zstd::libzstd_static)
# Different between building with installed vs source boost
if(NOT TARGET Boost::headers)
  set(SKETCHER_DEPS ${SKETCHER_DEPS} Boost::assign Boost::beast Boost::bimap)
endif()
# Different between building with installed vs source zstd
if(NOT TARGET zstd::libzstd_static)
  add_library(zstd::libzstd_static ALIAS libzstd_static)
endif()

# Sketcher Library
set(LIB_TARGET ${CMAKE_PROJECT_NAME})
file(GLOB_RECURSE SOURCE_LIST CONFIGURE_DEPENDS ${SKETCHER_INC}/*.h
     ${SKETCHER_SRC}/*.h ${SKETCHER_SRC}/*.cpp)
set(CMAKE_AUTOUIC_SEARCH_PATHS ${SKETCHER_SRC}/ui)
add_library(${LIB_TARGET} STATIC ${SOURCE_LIST})
target_compile_options(${LIB_TARGET} PRIVATE ${SKETCHER_COMPILE_FLAGS})
target_include_directories(${LIB_TARGET} PRIVATE ${SKETCHER_LIB_INCLUDES})
target_compile_definitions(${LIB_TARGET} PRIVATE SKETCHER_STATIC_DEFINE)
target_link_libraries(${LIB_TARGET} PRIVATE ${SKETCHER_DEPS})

# Sketcher Executable
set(EXE_TARGET sketcher_app)
add_executable(${EXE_TARGET} src/app/main.cpp)
target_compile_options(${EXE_TARGET} PRIVATE ${SKETCHER_COMPILE_FLAGS})
target_include_directories(${EXE_TARGET} PRIVATE ${SKETCHER_LIB_INCLUDES})
target_compile_definitions(${EXE_TARGET} PRIVATE SKETCHER_STATIC_DEFINE)
target_link_libraries(${EXE_TARGET} PRIVATE ${LIB_TARGET})

# Sketcher Tests
enable_testing()

file(GLOB_RECURSE TEST_SOURCE_LIST CONFIGURE_DEPENDS ${SKETCHER_TEST}/*.cpp)
set(TEST_DEP_LIBS ${LIB_TARGET} Boost::unit_test_framework ${SKETCHER_DEPS})
foreach(TEST_SOURCE ${TEST_SOURCE_LIST})
  cmake_path(GET TEST_SOURCE STEM TEST_TARGET)
  add_executable(${TEST_TARGET} ${TEST_SOURCE})
  target_compile_options(${TEST_TARGET} PRIVATE ${SKETCHER_COMPILE_FLAGS})
  target_include_directories(${TEST_TARGET} PRIVATE ${SKETCHER_LIB_INCLUDES})
  target_compile_definitions(${TEST_TARGET} PRIVATE SKETCHER_STATIC_DEFINE)
  target_link_libraries(${TEST_TARGET} PRIVATE ${TEST_DEP_LIBS})
  add_test(NAME ${TEST_TARGET} COMMAND ${TEST_TARGET})
  # Mirror the source directory for test collection
  cmake_path(RELATIVE_PATH TEST_SOURCE OUTPUT_VARIABLE TEST_PATH)
  cmake_path(GET TEST_PATH PARENT_PATH PREFIX_PATH)
  set_target_properties(${TEST_TARGET} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                  ${PREFIX_PATH})
endforeach()
