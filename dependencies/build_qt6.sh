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
install_dir=$(pwd)/dependencies/${platform}_${sha}_qt6

version_prefix=$(echo ${version} | cut -d. -f1-2)
num_cores=$(python3 -c 'import multiprocessing; print(multiprocessing.cpu_count())')

opts="-prefix ${install_dir} -release -static -nomake examples -nomake tests -make tools -no-opengl"

rm -fr qt*everywhere-src* || true

for lib in base svg; do
    if [ ! -f qt6-${lib}-${version}.tar.xz ]; then
        curl -L -o qt6-${lib}-${version}.tar.xz https://download.qt.io/archive/qt/${version_prefix}/${version}/submodules/qt${lib}-everywhere-opensource-src-${version}.tar.xz
    fi
    tar xf qt6-${lib}-${version}.tar.xz
    pushd qt${lib}-everywhere-src-${version}
    cmake \
        -Wno-dev \
        -G Ninja \
        -DQt6_DIR=${install_dir}/lib/cmake/Qt6 \
        -DBUILD_SHARED_LIBS=OFF \
        -DCMAKE_INSTALL_PREFIX=${install_dir} \
        -DQT_BUILD_EXAMPLES=FALSE \
        -DQT_BUILD_TESTS=FALSE \
        -DQT_FEATURE_brotli=OFF \
        -DQT_FEATURE_sql=OFF \
        -DCMAKE_BUILD_TYPE=Release \
        -DINPUT_opengl=no \
        $(pwd)
    cmake --build . -j${num_cores} && cmake --install .
    popd
done





