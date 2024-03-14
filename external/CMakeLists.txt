cmake_minimum_required(VERSION 3.24)
include(ExternalProject)

project(schrodinger_sketcher_external)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MSVC)
  # Workaround until boost 1.83 which fixes the issue with std::unary_function
  add_definitions(-D_HAS_AUTO_PTR_ETC=0)
endif()

# Collect dependency versions
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/versions.json VERSIONS)
string(JSON BOOST_VERSION GET ${VERSIONS} boost)
string(JSON FMT_VERSION GET ${VERSIONS} fmt)
string(JSON QT_VERSION GET ${VERSIONS} qt)
string(JSON RDKIT_VERSION GET ${VERSIONS} rdkit)
string(JSON ZSTD_VERSION GET ${VERSIONS} zstd)

ExternalProject_Add(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG ${FMT_VERSION}
  GIT_SHALLOW TRUE
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/fmt-${FMT_VERSION}
    -DCMAKE_BUILD_TYPE=Release -DFMT_DOC=OFF -DFMT_TEST=OFF -DFMT_INSTALL=ON)

ExternalProject_Add(
  zstd
  GIT_REPOSITORY https://github.com/facebook/zstd.git
  GIT_TAG v${ZSTD_VERSION}
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR build/cmake
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/zstd-${ZSTD_VERSION}
    -DCMAKE_BUILD_TYPE=Release -DZSTD_BUILD_SHARED=OFF
    -DZSTD_BUILD_PROGRAMS=OFF)

set(BOOST_INCLUDE_LIBS
    # RDKit required libraries
    assign
    crc
    format
    flyweight
    graph
    iostreams
    multi_array
    multiprecision
    program_options
    serialization
    # rdkit_extensions required libraries
    beast
    bimap
    filesystem
    test)
# Use generator expression to sidestep CMAKE_ARGS removing semicolons
string(REPLACE ";" "$<SEMICOLON>" BOOST_INCLUDE_LIBS "${BOOST_INCLUDE_LIBS}")
ExternalProject_Add(
  boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG boost-${BOOST_VERSION}
  GIT_SHALLOW TRUE
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/boost-${BOOST_VERSION}
    -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF
    -DBOOST_INCLUDE_LIBRARIES=${BOOST_INCLUDE_LIBS})

string(REPLACE "." "_" RDKIT_GIT_TAG ${RDKIT_VERSION})
ExternalProject_Add(
  rdkit
  GIT_REPOSITORY https://github.com/rdkit/rdkit.git
  GIT_TAG Release_${RDKIT_GIT_TAG}
  GIT_SHALLOW TRUE
  PATCH_COMMAND git apply ${CMAKE_CURRENT_SOURCE_DIR}/rdkit.patch
  UPDATE_DISCONNECTED 1 # required to not apply patch again once populated
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/rdkit-${RDKIT_VERSION}
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}/boost-${BOOST_VERSION}
    -DBUILD_SHARED_LIBS=OFF
    -DBoost_USE_STATIC_LIBS=ON
    -DRDK_BUILD_CPP_TESTS=OFF
    -DRDK_BUILD_DESCRIPTORS3D=OFF
    -DRDK_BUILD_FREETYPE_SUPPORT=OFF
    -DRDK_BUILD_INCHI_SUPPORT=ON
    -DRDK_BUILD_MAEPARSER_SUPPORT=ON
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF
    -DRDK_BUILD_SLN_SUPPORT=OFF
    -DRDK_BUILD_XYZ2MOL_SUPPORT=ON
    -DRDK_INSTALL_COMIC_FONTS=OFF
    -DRDK_INSTALL_INTREE=OFF
    -DRDK_INSTALL_STATIC_LIBS=OFF
    -DRDK_OPTIMIZE_POPCNT=OFF
    -DRDK_USE_BOOST_STACKTRACE=OFF)

ExternalProject_Add(
  qt
  GIT_REPOSITORY https://github.com/qt/qt5.git
  GIT_TAG v${QT_VERSION}-lts-lgpl
  GIT_SUBMODULES qtbase qtsvg
  GIT_SHALLOW TRUE
  PATCH_COMMAND cd qtsvg && git apply ${CMAKE_CURRENT_SOURCE_DIR}/qtsvg.patch
  UPDATE_DISCONNECTED 1 # required to not apply patch again once populated
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/qt-${QT_VERSION}
             -DCMAKE_BUILD_TYPE=Release
             -DBUILD_SHARED_LIBS=OFF
             -DQT_BUILD_EXAMPLES=OFF
             -DQT_BUILD_SUBMODULES=qtbase$<SEMICOLON>qtsvg
             -DQT_BUILD_TESTS=OFF
             -DQT_FEATURE_brotli=OFF
             -DQT_FEATURE_sql=OFF
             -DINPUT_opengl=no)
