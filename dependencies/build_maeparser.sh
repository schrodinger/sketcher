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
install_dir=$(pwd)/dependencies/${platform}_${sha}_maeparser
boost_root=$(pwd)/dependencies/${platform}_${sha}_boost
zlib_root=$(pwd)/dependencies/${platform}_${sha}_zlib

if [ ! -f maeparser-${version}.tar.gz ]; then
    curl -L -o maeparser-${version}.tar.gz https://github.com/schrodinger/maeparser/archive/refs/tags/v${version}.tar.gz
fi

rm -fr maeparser-${version} || true

tar xf maeparser-${version}.tar.gz && cd maeparser-${version}

cmake \
    -Wno-dev \
    -G Ninja \
    -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${install_dir} \
    -DMAEPARSER_BUILD_TESTS=OFF \
    -DMAEPARSER_BUILD_SHARED_LIBS=OFF \
    -DBOOST_ROOT=${boost_root} \
    -DBoost_USE_DEBUG_RUNTIME=OFF \
    -DBoost_USE_MULTITHREADED=OFF \
    -DBoost_USE_STATIC_RUNTIME=ON \
    -DZLIB_ROOT=${zlib_root} \
    .

cd build && ninja install



