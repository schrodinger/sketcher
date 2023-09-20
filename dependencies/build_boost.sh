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
boost_build_opts="address-model=64 architecture=x86 link=static runtime-link=static threading=single variant=release"
num_cores=$(python3 -c 'import multiprocessing; print(multiprocessing.cpu_count())')

if [ "${platform}" == "win32" ]; then
    msvc_toolset_major=${VC_REDIST_VERSION:0:2}
    msvc_toolset_minor=${VC_REDIST_VERSION:2:1}
    boost_build_opts="${boost_build_opts} toolset=msvc-${msvc_toolset_major}.${msvc_toolset_minor}"
else
    boost_build_opts="${boost_build_opts} cxxflags=-fPIC cflags=-fPIC"
fi

if [ ! -f boost-${boost_version}.tar.gz ]; then
    curl -L -o boost-${boost_version}.tar.gz https://boostorg.jfrog.io/artifactory/main/release/${version}/source/boost_${boost_version}.tar.gz
fi

rm -fr boost_${boost_version} || true

tar xf boost-${boost_version}.tar.gz && cd boost_${boost_version}

./bootstrap.sh --with-libraries=${boost_libs} --prefix=${install_dir}

./b2 --prefix=${install_dir} ${boost_build_opts} -s NO_ZSTD=1 -s NO_BZIP2=1 -j ${num_cores} install