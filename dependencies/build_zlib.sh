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
install_dir=$(pwd)/dependencies/${platform}_${sha}_zlib

if [ ! -f zlib-${version}.tar.gz ]; then
    curl -L -o zlib-${version}.tar.gz https://github.com/madler/zlib/archive/refs/tags/v${version}.tar.gz
fi

rm -fr zlib-${version} || true

tar xf zlib-${version}.tar.gz && cd zlib-${version}

# zlib shared lib build cannot be disabled.
# Let's just not install it

if [ "${platform}" == "darwin" ]; then
    sed_arg='""'
fi

sed -i ${sed_arg:-} 's/install(TARGETS zlib zlibstatic/install(TARGETS zlibstatic/' CMakeLists.txt

cmake \
    -Wno-dev \
    -G Ninja \
    -B build \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX=${install_dir} \
    -DBUILD_SHARED_LIBS=OFF \
    .

cd build && ninja install

