#!/bin/bash

set -x -e -u -o pipefail

if [ $# -lt 3 ]; then
    cmd=$(basename $0)
    echo "Usage: ${cmd} <platform> <version> <commit sha>"
    echo "Actual call was: ${cmd} $*"
    exit 1
fi

platform=$1
version=$2
sha=$3
install_dir=$(pwd)/dependencies/${platform}_${sha}_rdkit
boost_root=$(pwd)/dependencies/${platform}_${sha}_boost
maeparser_root=$(pwd)/dependencies/${platform}_${sha}_maeparser
zlib_root=$(pwd)/dependencies/${platform}_${sha}_zlib

patch_file=$(pwd)/dependencies/rdk_static_libs_only.patch

if [ ! -f rdkit-${version}.tar.gz ]; then
    curl -L -o rdkit-${version}.tar.gz https://github.com/rdkit/rdkit/archive/${version}.tar.gz
fi

rm -fr rdkit-${version} || true

tar xf rdkit-${version}.tar.gz && cd rdkit-${version}

# Patch to avoid building the shared libraries. This helps reducing
# Besidese the patch, RDK_INSTALL_STATIC_LIBS=OFF should be set.
patch -p1 < ${patch_file}

# Required to find pre-built maeparser
export CMAKE_PREFIX_PATH=${maeparser_root}

cmake \
    -Wno-dev \
    -G Ninja \
    -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${install_dir} \
    -DBOOST_ROOT=${boost_root} \
    -DZLIB_ROOT=${zlib_root} \
    -Dmaeparser_DIR=${maeparser_root}/lib/cmake \
    -DBoost_USE_DEBUG_RUNTIME=OFF \
    -DBoost_USE_MULTITHREADED=OFF \
    -DBoost_USE_STATIC_LIBS=ON \
    -DBoost_USE_STATIC_RUNTIME=ON \
    -DRDK_INSTALL_STATIC_LIBS=OFF \
    -DRDK_BUILD_CPP_TESTS=OFF \
    -DRDK_BUILD_DESCRIPTORS3D=OFF \
    -DRDK_BUILD_FREETYPE_SUPPORT=OFF \
    -DRDK_BUILD_INCHI_SUPPORT=ON \
    -DRDK_BUILD_MAEPARSER_SUPPORT=ON \
    -DRDK_BUILD_PYTHON_WRAPPERS=OFF \
    -DRDK_BUILD_SLN_SUPPORT=OFF \
    -DRDK_BUILD_XYZ2MOL_SUPPORT=ON \
    -DRDK_INSTALL_COMIC_FONTS=OFF \
    -DRDK_INSTALL_INTREE=OFF \
    -DRDK_USE_BOOST_STACKTRACE=OFF \
    .

cd build && ninja install
