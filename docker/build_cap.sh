#!/usr/bin/env bash

cd ${PREFIX}/source && \
git clone https://github.com/dalg24/cap-data.git && \
mkdir -p ${PREFIX}/build/cap && \
cd ${PREFIX}/build/cap && \
cmake \
    -D CMAKE_INSTALL_PREFIX=/opt/cap \
    -D BUILD_SHARED_LIB=ON \
    -D MPI_INSTALL_DIR=${MPI_DIR} \
    -D PYTHON_INSTALL_DIR=${PYTHON_DIR} \
    -D BOOST_INSTALL_DIR=${BOOST_DIR} \
    -D DEAL_II_INSTALL_DIR=${DEAL_II_DIR} \
    -D CAP_DATA_DIR=${PREFIX}/source/cap-data \
    ${PREFIX}/source/cap && \
make install && \
cd python && \
ctest -V
