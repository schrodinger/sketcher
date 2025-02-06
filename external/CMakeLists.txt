cmake_minimum_required(VERSION 3.24)
include(ExternalProject)

project(schrodinger_sketcher_external)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Collect dependency versions
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/versions.json VERSIONS)
string(JSON BOOST_VERSION GET ${VERSIONS} boost)
string(JSON FMT_VERSION GET ${VERSIONS} fmt)
string(JSON QT_VERSION GET ${VERSIONS} qt)
string(JSON RDKIT_VERSION GET ${VERSIONS} rdkit)
string(JSON ZLIB_VERSION GET ${VERSIONS} zlib)
string(JSON ZSTD_VERSION GET ${VERSIONS} zstd)

# Ensure args passed to ExternalProject_Add's CMAKE_ARGS support semicolons
function(set_cmake_arg VAR_NAME)
  set(CMAKE_ARG ${ARGN})
  string(REPLACE ";" "$<SEMICOLON>" CMAKE_ARG "${CMAKE_ARG}")
  set(${VAR_NAME}
      "${CMAKE_ARG}"
      PARENT_SCOPE)
endfunction()

if(LINUX)
    set(EXTRA_QT_CMAKE_ARGS
        -DFEATURE_xlib=ON
        -DFEATURE_xcb=ON
        -DFEATURE_xcb_xlib=ON
        -DFEATURE_xkbcommon=ON
    )
else()
    set(EXTRA_QT_CMAKE_ARGS "")
endif()

ExternalProject_Add(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG ${FMT_VERSION}
  GIT_SHALLOW TRUE
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/fmt-${FMT_VERSION}
    -DCMAKE_BUILD_TYPE=Release -DFMT_DOC=OFF -DFMT_TEST=OFF -DFMT_INSTALL=ON)

ExternalProject_Add(
  zlib
  GIT_REPOSITORY https://github.com/madler/zlib.git
  GIT_TAG v${ZLIB_VERSION}
  GIT_SHALLOW TRUE
  PATCH_COMMAND git apply ${CMAKE_CURRENT_SOURCE_DIR}/zlib.patch
  UPDATE_DISCONNECTED 1 # required to not apply patch again once populated
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/zlib-${ZLIB_VERSION}
    -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF)

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

set_cmake_arg(
  BOOST_INCLUDE_LIBS
  # RDKit required libraries
  assign
  crc
  format
  flyweight
  graph
  iostreams
  json
  multi_array
  multiprecision
  program_options
  serialization
  # rdkit_extensions required libraries
  beast
  bimap
  filesystem
  test)
set_cmake_arg(
  BOOST_CMAKE_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/zlib-${ZLIB_VERSION}
  ${CMAKE_CURRENT_BINARY_DIR}/zstd-${ZSTD_VERSION})
ExternalProject_Add(
  boost
  GIT_REPOSITORY https://github.com/boostorg/boost.git
  GIT_TAG boost-${BOOST_VERSION}
  GIT_SHALLOW TRUE
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/boost-${BOOST_VERSION}
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_PREFIX_PATH=${BOOST_CMAKE_PREFIX}
    -DBUILD_SHARED_LIBS=OFF
    -DBOOST_INCLUDE_LIBRARIES=${BOOST_INCLUDE_LIBS}
    -DBOOST_RUNTIME_LINK=static
    -DBOOST_IOSTREAMS_ENABLE_BZIP2=OFF
    -DBOOST_IOSTREAMS_ENABLE_LZMA=OFF)

string(REPLACE "." "_" RDKIT_GIT_TAG ${RDKIT_VERSION})
set_cmake_arg(
  RDKIT_CMAKE_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/boost-${BOOST_VERSION}
  ${CMAKE_CURRENT_BINARY_DIR}/zlib-${ZLIB_VERSION})
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
    -DCMAKE_PREFIX_PATH=${RDKIT_CMAKE_PREFIX}
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

set_cmake_arg(QT_SUBMODULES qtbase qtsvg)
ExternalProject_Add(
  qt
  GIT_REPOSITORY https://github.com/qt/qt5.git
  GIT_TAG v${QT_VERSION}
  GIT_SUBMODULES qtbase qtsvg
  GIT_SHALLOW TRUE
  CMAKE_ARGS
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/qt-${QT_VERSION}
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}/zlib-${ZLIB_VERSION}
    -DCMAKE_POLICY_DEFAULT_CMP0177:STRING=NEW
    -DCMAKE_POLICY_DEFAULT_CMP0174:STRING=NEW
    -DBUILD_SHARED_LIBS=OFF
    -DQT_BUILD_EXAMPLES=OFF
    -DQT_BUILD_SUBMODULES=${QT_SUBMODULES}
    -DQT_BUILD_TESTS=OFF
    -DQT_FEATURE_brotli=OFF
    -DQT_FEATURE_sql=OFF
    -DINPUT_opengl=no
    ${EXTRA_QT_CMAKE_ARGS})
