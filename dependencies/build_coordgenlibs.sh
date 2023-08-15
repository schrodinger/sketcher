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
install_dir=$(pwd)/dependencies/${platform}_${sha}_coordgenlibs

if [ ! -f coordgenlibs-${version}.tar.gz ]; then
    curl -L -o coordgenlibs-${version}.tar.gz https://github.com/schrodinger/coordgenlibs/archive/refs/tags/v${version}.tar.gz
fi

rm -fr coordgenlibs-${version} || true

tar xf coordgenlibs-${version}.tar.gz && cd coordgenlibs-${version}

cmake \
    -Wno-dev \
    -G Ninja \
    -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_INSTALL_PREFIX=${install_dir} \
    -DCOORDGEN_BUILD_EXAMPLE=OFF \
    -DCOORDGEN_BUILD_TESTS=OFF \
    -DCOORDGEN_BUILD_SHARED_LIBS=OFF \
    -DCOORDGEN_USE_MAEPARSER=OFF \
    .

cd build && ninja install



