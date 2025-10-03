#!/bin/bash

set -euo pipefail

export INSTALL_DIR=$1

if [[ -z "$INSTALL_DIR" ]]; then
    echo "Must give filename for the tarball" >&2
    exit 1
fi
rm -fr ${INSTALL_DIR} || true
mkdir -p ${INSTALL_DIR}

echo -e "Building WASM distribution in ${INSTALL_DIR}"
echo -e "\tRunning in:" `pwd` 

ls -al ./build
cp -v ./build/sketcher_app.js ${INSTALL_DIR}
cp -v ./build/sketcher_app.wasm ${INSTALL_DIR}
cp -v ${QT_DIR}/plugins/platforms/qtloader.js ${INSTALL_DIR}
cp -vr ./wasm/public/. ${INSTALL_DIR}

export md5Wasm=$(md5sum ./build/sketcher_app.wasm | awk '{ print $1 }')
perl -i -pe 's/sketcher_app\.wasm/sketcher_app.wasm?cache_bust=$ENV{md5Wasm}/g' ${INSTALL_DIR}/sketcher_app.js

export md5Qt=$(md5sum ${INSTALL_DIR}/qtloader.js | awk '{ print $1 }')
export md5Js=$(md5sum ./build/sketcher_app.js | awk '{ print $1 }')
perl -i -pe 's/qtloader\.js/qtloader.js?cache_bust=$ENV{md5Qt}/g' ${INSTALL_DIR}/wasm_shell.html
perl -i -pe 's/sketcher_app\.js/sketcher_app.js?cache_bust=$ENV{md5Js}/g' ${INSTALL_DIR}/wasm_shell.html

echo -e "Built WASM distribution in ${INSTALL_DIR}\n"

echo -e "Creating ${INSTALL_DIR}.tar.gz"
rm -f ${INSTALL_DIR}.tar.gz || true
tar czvf ${INSTALL_DIR}.tar.gz  ${INSTALL_DIR}
echo -e "Successfully created ${INSTALL_DIR}.tar.gz"