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
install_dir=$(pwd)/dependencies/${platform}_${sha}_fmt

if [ ! -f fmt-${version}.tar.gz ]; then
    curl -L -o fmt-${version}.tar.gz https://github.com/fmtlib/fmt/archive/refs/tags/${version}.tar.gz
fi

rm -fr fmt-${version} || true

tar xf fmt-${version}.tar.gz && cd fmt-${version}

cmake \
    -Wno-dev \
    -G Ninja \
    -B build \
    -DBUILD_SHARED_LIBS=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${install_dir} \
    -DFMT_DOC=OFF \
    -DFMT_TEST=OFF \
    -DFMT_INSTALL=ON \
    .

cd build && ninja install

