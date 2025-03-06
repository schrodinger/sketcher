#!/bin/bash

export INSTALL_DIR=$1

rm -fr ${INSTALL_DIR} || true
mkdir -p ${INSTALL_DIR}

cp ./build/sketcher.js ${INSTALL_DIR}
cp ./build/sketcher.wasm ${INSTALL_DIR}
cp ${QT_DIR}/plugins/platforms/qtloader.js ${INSTALL_DIR}
cp -r ./wasm/public/. ${INSTALL_DIR}

export md5Wasm=$(md5sum sketcher.wasm | awk '{ print $1 }')
export md5Js=$(md5sum sketcher.js | awk '{ print $1 }')

perl -i -pe 's/\.wasm/.wasm?cache_bust=$ENV{md5Wasm}/g' ${INSTALL_DIR}/qtloader.js
perl -i -pe 's/"\.js"/".js?cache_bust=$ENV{md5Js}"/g' ${INSTALL_DIR}/qtloader.js

export md5Qt=$(md5sum ${INSTALL_DIR}/qtloader.js | awk '{ print $1 }')
perl -i -pe 's/qtloader\.js/qtloader.js?cache_bust=$ENV{md5Qt}/g' ${INSTALL_DIR}/wasm_shell.html

echo -e "\n\nBuilt WASM distibution in ${INSTALL_DIR}\n"

rm -f ${INSTALL_DIR}.tar.gz || true
tar czvf ${INSTALL_DIR}.tar.gz  ${INSTALL_DIR}
echo -e "\n\nSuccessfully created ${INSTALL_DIR}.tar.gz\n"