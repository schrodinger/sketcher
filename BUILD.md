## How to build the Schrödinger 2D Sketcher

If you are a Schrödinger developer using the Core Suite build environment,
simply use `buildinger sketcher`.

To build outside of the Schrödinger build environment, the recommended
procedure on Linux is to create an Ubuntu 22.04 Docker image that includes the
following packages:

    build-essential
    curl
    git
    locales
    pkg-config
    python3
    python3-dev
    python3-pip
    python3-venv
    sudo
    vim
    wget
    libfontconfig1-dev
    libfreetype-dev
    libx11-dev
    libx11-xcb-dev
    libxcb-cursor-dev
    libxcb-glx0-dev
    libxcb-icccm4-dev
    libxcb-image0-dev
    libxcb-keysyms1-dev
    libxcb-randr0-dev
    libxcb-render-util0-dev
    libxcb-shape0-dev
    libxcb-shm0-dev
    libxcb-sync-dev
    libxcb-util-dev
    libxcb-xfixes0-dev
    libxcb-xinerama0-dev
    libxcb-xkb-dev
    libxcb1-dev
    libxext-dev
    libxfixes-dev
    libxi-dev
    libxkbcommon-dev
    libxkbcommon-x11-dev
    libxrender-dev
    xvfb

Once in the container and in the sketcher directory, create a virtualenv to get
the correct versions of cmake and ninja (and pytest, for the unit tests):

    python3 -m venv .venv
    . .venv/bin/activate
    pip install cmake==3.24 ninja==1.13.0 pytest pytest-cpp pytest-xdist

Next build the external dependencies:

    cmake -B ext_bld -S external -G Ninja -DCMAKE_BUILD_TYPE=Release

    # The first five are quick to build and independent of each other, so
    # they can be built together.
    cmake --build ext_bld --target fmt zlib zstd eigen sqlite
    cmake --build ext_bld --target boost
    cmake --build ext_bld --target rdkit
    cmake --build ext_bld --target qt

    CMAKE_PREFIX_PATH=`echo ext_bld/*-[0-9]* | sed 's/ /;/g'`

And finally build the sketcher itself:

    cmake -B build -S . -G Ninja -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_FIND_PACKAGE_PREFER_CONFIG=ON \
      -DCMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH

    cmake --build build --config Release

To install files to their final locations:

    cmake --install build --config Release --prefix build

To run the unit tests:

    cd build
    SKETCHER_SOURCE_DIR=$PWD/.. xvfb-run -a pytest -n auto
