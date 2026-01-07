cmake_minimum_required(VERSION 3.24)
include(ExternalProject)

project(schrodinger_sketcher_external)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Collect dependency versions
file(READ ${CMAKE_CURRENT_SOURCE_DIR}/versions.json VERSIONS)
string(JSON BOOST_VERSION GET ${VERSIONS} boost)
string(JSON EIGEN_VERSION GET ${VERSIONS} eigen)
string(JSON FMT_VERSION GET ${VERSIONS} fmt)
string(JSON QT_VERSION GET ${VERSIONS} qt)
string(JSON RDKIT_VERSION GET ${VERSIONS} rdkit)
string(JSON SQLITE_VERSION GET ${VERSIONS} sqlite)
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

# Generate URLS for downloading the SQLite amalgamation, which contains all of
# the SQLite source code concatenated into a single header and a single
# implementation file. The URL contains the release year, but there's no trivial
# way to determine the year from the version number, so we just try a range of
# years until we find the correct one. This could fail for edge cases (e.g. the
# filenames for versions 3.5.0, 3.50.0, and 3.500.0 will all be the same, so we
# might find the wrong one if our range of years is way too large), but it
# should be good enough for our needs.
function(generate_sqlite_urls VAR_NAME)
  # pull out the digits from the SQLite version
  string(REGEX REPLACE "[^0-9]" "" SQLITE_VERS_DIGITS "${SQLITE_VERSION}")
  # the file name is padded with zeroes to seven digits
  string(LENGTH "${SQLITE_VERS_DIGITS}" SQLITE_VERS_NUM_DIGITS)
  math(EXPR NUM_ZEROS "7 - ${SQLITE_VERS_NUM_DIGITS}")
  string(REPEAT "0" ${NUM_ZEROS} TRAILING_ZEROS)
  foreach(YEAR RANGE 2023 2035)
    list(
      APPEND
      URLS_FOR_YEARS
      "https://www.sqlite.org/${YEAR}/sqlite-amalgamation-${SQLITE_VERS_DIGITS}${TRAILING_ZEROS}.zip"
    )
  endforeach(YEAR)
  list(JOIN URLS_FOR_YEARS ";" JOINED_URLS)
  set(${VAR_NAME}
      "${JOINED_URLS}"
      PARENT_SCOPE)
endfunction()

if(EMSCRIPTEN)
  set(EXTRA_BOOST_CMAKE_ARGS
      -DCMAKE_CXX_FLAGS=-fdeclspec\ -pthread\ -fexceptions\ -sDISABLE_EXCEPTION_CATCHING=0
  )
  set(EXTRA_QT_CMAKE_ARGS
      -DQT_HOST_PATH=${CMAKE_CURRENT_BINARY_DIR}/host/qt-${QT_VERSION})
  set(EXTRA_RDKIT_CMAKE_ARGS
      -DCMAKE_CXX_FLAGS=-fdeclspec\ -pthread\ -fexceptions\ -sDISABLE_EXCEPTION_CATCHING=0
      -DCMAKE_FIND_ROOT_PATH=${CMAKE_CURRENT_BINARY_DIR})
else()
  set(EXTRA_QT_CMAKE_ARGS -DINPUT_opengl=no)
  if(LINUX)
    list(APPEND EXTRA_QT_CMAKE_ARGS -DFEATURE_xlib=ON -DFEATURE_xcb=ON
         -DFEATURE_xcb_xlib=ON -DFEATURE_xkbcommon=ON)
  endif()
endif()

set(COMMON_EXTERNAL_PROJECT_CMAKE_ARGS
    -DCMAKE_CXX_STANDARD=${CMAKE_CXX_STANDARD} -DCMAKE_CXX_STANDARD_REQUIRED=ON
    -DCMAKE_BUILD_TYPE=Release)
if(CMAKE_TOOLCHAIN_FILE)
  # Emscripten builds specify a toolchain file and a cross-compiling emulator
  list(APPEND COMMON_EXTERNAL_PROJECT_CMAKE_ARGS
       -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
       -DCMAKE_CROSSCOMPILING_EMULATOR=${CMAKE_CROSSCOMPILING_EMULATOR})
endif()
if(CMAKE_FIND_ROOT_PATH)
  # make sure we can always find the dependencies that we've built, even if
  # emscripten is trying to prevent us from finding system libraries
  list(APPEND COMMON_EXTERNAL_PROJECT_CMAKE_ARGS
       -DCMAKE_FIND_ROOT_PATH=${CMAKE_CURRENT_BINARY_DIR})
endif()
if(MSVC)
  # force MSVC to use the static linking runtime library instead of the default
  # shared version
  list(APPEND COMMON_EXTERNAL_PROJECT_CMAKE_ARGS
       -DCMAKE_MSVC_RUNTIME_LIBRARY=MultiThreaded$<$<CONFIG:Debug>:Debug>)
endif()

# Build fmt library
ExternalProject_Add(
  fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG ${FMT_VERSION}
  GIT_SHALLOW TRUE
  CMAKE_ARGS
    ${COMMON_EXTERNAL_PROJECT_CMAKE_ARGS}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/fmt-${FMT_VERSION}
    -DFMT_DOC=OFF -DFMT_TEST=OFF -DFMT_INSTALL=ON)

# Build zlib library
ExternalProject_Add(
  zlib
  GIT_REPOSITORY https://github.com/madler/zlib.git
  GIT_TAG v${ZLIB_VERSION}
  GIT_SHALLOW TRUE
  PATCH_COMMAND git reset --hard HEAD && git apply
                ${CMAKE_CURRENT_SOURCE_DIR}/zlib_no_shared.patch
  CMAKE_ARGS
    ${COMMON_EXTERNAL_PROJECT_CMAKE_ARGS}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/zlib-${ZLIB_VERSION}
    -DBUILD_SHARED_LIBS=OFF -DZLIB_BUILD_EXAMPLES=OFF)

# Build zstd library
ExternalProject_Add(
  zstd
  GIT_REPOSITORY https://github.com/facebook/zstd.git
  GIT_TAG v${ZSTD_VERSION}
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR build/cmake
  CMAKE_ARGS
    ${COMMON_EXTERNAL_PROJECT_CMAKE_ARGS}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/zstd-${ZSTD_VERSION}
    -DZSTD_BUILD_SHARED=OFF -DZSTD_BUILD_PROGRAMS=OFF)

# Build Eigen library
ExternalProject_Add(
  eigen
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG ${EIGEN_VERSION}
  GIT_SHALLOW TRUE
  CMAKE_ARGS
    ${COMMON_EXTERNAL_PROJECT_CMAKE_ARGS}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/eigen-${EIGEN_VERSION})

# Build SQLite library
generate_sqlite_urls(SQLITE_URLS)
ExternalProject_Add(
  sqlite
  URL ${SQLITE_URLS}
  # SQLite doesn't work with CMake out of the box, so copy in our CMake files
  PATCH_COMMAND ${CMAKE_COMMAND} -E copy_directory
                ${CMAKE_CURRENT_SOURCE_DIR}/sqlite_cmake <SOURCE_DIR>
  CMAKE_ARGS
    ${COMMON_EXTERNAL_PROJECT_CMAKE_ARGS}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/sqlite-${SQLITE_VERSION}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}/zlib-${ZLIB_VERSION})

# Build Boost library
set_cmake_arg(
  BOOST_INCLUDE_LIBS
  # RDKit required libraries
  assign
  crc
  format
  flyweight
  graph
  geometry
  iostreams
  json
  locale
  multi_array
  multiprecision
  program_options
  rational
  serialization
  # rdkit_extensions required libraries
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
    ${COMMON_EXTERNAL_PROJECT_CMAKE_ARGS}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/boost-${BOOST_VERSION}
    -DCMAKE_PREFIX_PATH=${BOOST_CMAKE_PREFIX}
    -DBUILD_SHARED_LIBS=OFF
    -DBOOST_INCLUDE_LIBRARIES=${BOOST_INCLUDE_LIBS}
    -DBOOST_RUNTIME_LINK=static
    -DBOOST_IOSTREAMS_ENABLE_BZIP2=OFF
    -DBOOST_IOSTREAMS_ENABLE_LZMA=OFF
    -DBOOST_IOSTREAMS_ENABLE_ZSTD=ON
    -DBOOST_IOSTREAMS_ZSTD_TARGET=zstd::libzstd_static
    ${EXTRA_BOOST_CMAKE_ARGS})

# building boost::beast will fail on wasm builds because it pulls in
# boost::context (which contains assembly), but we can install beast as a
# header-only library and use it that way
set(BOOST_BEAST_HEADERS_PATH <SOURCE_DIR>/libs/beast/include)
set(BOOST_INCLUDE_DIR ${CMAKE_CURRENT_BINARY_DIR}/boost-${BOOST_VERSION})
if(WIN32)
  string(APPEND BOOST_BEAST_HEADERS_PATH "/*")
  # we need to append /include/boost-1_81 directory to the include directory
  # path, but BOOST_VERSION uses periods, not underscores
  string(REPLACE "." ";" BOOST_VERSION_PARTS ${BOOST_VERSION})
  list(POP_BACK BOOST_VERSION_PARTS)
  list(JOIN BOOST_VERSION_PARTS "_" BOOST_VERSION_USING_UNDERSCORE)
  string(APPEND BOOST_INCLUDE_DIR
         "/include/boost-${BOOST_VERSION_USING_UNDERSCORE}")
else()

endif()
ExternalProject_Add_Step(
  boost install_beast_header_only
  COMMENT "Copying boost::beast headers"
  DEPENDEES install
  ALWAYS ON
  COMMAND cp -rv ${BOOST_BEAST_HEADERS_PATH} ${BOOST_INCLUDE_DIR})

# Build RDKit library
string(REPLACE "." "_" RDKIT_GIT_TAG ${RDKIT_VERSION})
set_cmake_arg(
  RDKIT_CMAKE_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/boost-${BOOST_VERSION}
  ${CMAKE_CURRENT_BINARY_DIR}/zlib-${ZLIB_VERSION}
  ${CMAKE_CURRENT_BINARY_DIR}/zstd-${ZSTD_VERSION})
set(RDKIT_PATCH_COMMAND git reset --hard HEAD)
# Prevent RDKit from trying to build shared libs
list(APPEND RDKIT_PATCH_COMMAND && git apply
     ${CMAKE_CURRENT_SOURCE_DIR}/rdkit_static_build.patch)
list(APPEND RDKIT_PATCH_COMMAND && git apply
     ${CMAKE_CURRENT_SOURCE_DIR}/rdkit_extern_patch.patch)
list(APPEND RDKIT_PATCH_COMMAND && git apply
     ${CMAKE_CURRENT_SOURCE_DIR}/rdkit_pr_9006_9016.patch)
if(WIN32)
  set(RDKIT_INSTALL_STATIC_LIBS_VAL ON)
else()
  set(RDKIT_INSTALL_STATIC_LIBS_VAL OFF)
endif()
ExternalProject_Add(
  rdkit
  GIT_REPOSITORY https://github.com/rdkit/rdkit.git
  GIT_TAG Release_${RDKIT_GIT_TAG}
  GIT_SHALLOW TRUE
  PATCH_COMMAND ${RDKIT_PATCH_COMMAND}
  CMAKE_ARGS
    ${COMMON_EXTERNAL_PROJECT_CMAKE_ARGS}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/rdkit-${RDKIT_VERSION}
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
    -DRDK_INSTALL_STATIC_LIBS=${RDKIT_INSTALL_STATIC_LIBS_VAL}
    -DRDK_OPTIMIZE_POPCNT=OFF
    -DRDK_USE_BOOST_STACKTRACE=OFF
    ${EXTRA_RDKIT_CMAKE_ARGS})

# Build Qt library
if(QT_HOST_BUILD)
  # For WASM builds, we also need a non-WASM build of Qt since some of the Qt
  # command line tools get run during the build process (e.g. uic, rcc, etc).
  # For those Qt builds, we don't need anything beyond qtbase.
  set(QT_GIT_SUBMODULES qtbase)
else()
  set(QT_GIT_SUBMODULES qtbase qtsvg)
endif()
set_cmake_arg(QT_BUILD_SUBMODULES ${QT_GIT_SUBMODULES})
ExternalProject_Add(
  qt
  GIT_REPOSITORY https://github.com/qt/qt5.git
  GIT_TAG v${QT_VERSION}
  GIT_SUBMODULES ${QT_GIT_SUBMODULES}
  GIT_SHALLOW TRUE
  CMAKE_ARGS
    ${COMMON_EXTERNAL_PROJECT_CMAKE_ARGS}
    -DCMAKE_INSTALL_PREFIX=${CMAKE_CURRENT_BINARY_DIR}/qt-${QT_VERSION}
    -DCMAKE_PREFIX_PATH=${CMAKE_CURRENT_BINARY_DIR}/zlib-${ZLIB_VERSION}
    -DCMAKE_POLICY_DEFAULT_CMP0177:STRING=NEW
    -DCMAKE_POLICY_DEFAULT_CMP0174:STRING=NEW
    -DBUILD_SHARED_LIBS=OFF
    -DQT_BUILD_EXAMPLES=OFF
    -DQT_BUILD_SUBMODULES=${QT_BUILD_SUBMODULES}
    -DQT_BUILD_TESTS=OFF
    -DQT_FEATURE_brotli=OFF
    -DQT_FEATURE_sql=OFF
    ${EXTRA_QT_CMAKE_ARGS})
