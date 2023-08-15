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
install_dir=$(pwd)/dependencies/${platform}_${sha}_boost

boost_version=${version//./_}
boost_libs="filesystem,iostreams,regex,serialization,system"
boost_build_opts="link=static runtime-link=static threading=single variant=release cxxflags=-fPIC cflags=-fPIC"
num_cores=$(python3 -c 'import multiprocessing; print(multiprocessing.cpu_count())')

if [ ! -f boost-${boost_version}.tar.gz ]; then
    curl -L -o boost-${boost_version}.tar.gz https://boostorg.jfrog.io/artifactory/main/release/${version}/source/boost_${boost_version}.tar.gz
fi

rm -fr boost_${boost_version} || true

tar xf boost-${boost_version}.tar.gz && cd boost_${boost_version}

./bootstrap.sh --with-libraries=${boost_libs} --prefix=${install_dir}

./b2 install ${boost_build_opts} -j ${num_cores}